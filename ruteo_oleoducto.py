"""
Ruteo óptimo de un oleoducto de recolección en Loma Campana (Vaca Muerta, Neuquén)
=================================================================================
Todo se descarga automáticamente de fuentes públicas (sin API keys):
  - DEM: Copernicus GLO-30 (AWS Open Data, 30 m)
  - Rutas, ductos, cursos de agua, zona urbana: OpenStreetMap (Overpass vía osmnx)

Instalación (Colab):  !pip install rasterio osmnx geopandas shapely scikit-image folium matplotlib
Uso:                  python ruteo_oleoducto.py            # datos reales
                      python ruteo_oleoducto.py --synthetic # prueba de lógica sin internet
"""
import json, sys, os
import numpy as np
from skimage.graph import route_through_array
from scipy.ndimage import distance_transform_edt

# ----------------------------- CONFIGURACIÓN -----------------------------
# A: coordenada de Loma Campana según Global Energy Monitor (WGS84)
# B: YPF Loma Campana - PTC (Planta de Tratamiento de Crudo), según Google Maps.
A = (-38.324, -68.831)   # (lat, lon)
B = (-38.355595, -68.712624)   # (lat, lon) YPF Loma Campana - PTC (Planta de Tratamiento de Crudo), Google Maps
BUFFER_M = 4000          # margen alrededor de A-B para el área de estudio
CRS_UTM = "EPSG:32719"   # WGS84 / UTM 19S
RES = 30                 # m por celda
NEAR_INFRA_M = 150       # distancia para considerar "paralelo" a infraestructura
ROAD_CROSS = 40          # costo extra por cruzar una ruta principal (tapada / perforación)
OUT = "outputs"

# ----------------------------- COSTOS -----------------------------------
def slope_percent(dem, res):
    gy, gx = np.gradient(dem, res)
    return np.hypot(gx, gy) * 100.0

def slope_cost(s):
    """Costo por pendiente (%): llano ~1, penalización creciente, muro >25%."""
    c = np.ones_like(s)
    c += (s / 5.0) ** 2                  # suave hasta ~10%
    c = np.where(s > 15, c + 60, c)      # obras de contención / zanjeo especial
    c = np.where(s > 25, c + 500, c)     # prácticamente prohibitivo (barda)
    return c

def build_cost(dem, res, infra_mask, water_mask, urban_mask, major_mask):
    s = slope_percent(dem, res)
    cost = slope_cost(s)
    d_infra = distance_transform_edt(~infra_mask) * res
    # Corredor: descuento solo si va PARALELO y desplazado (>= 1 celda), nunca sobre la calzada
    banda = (d_infra >= res) & (d_infra <= NEAR_INFRA_M)
    cost = np.where(banda, cost * 0.6, cost)
    cost = np.where(major_mask, cost + ROAD_CROSS, cost)         # cruce de ruta principal
    cost = np.where(water_mask, cost + 40, cost)   # cruce especial (finito: arroyos secos)
    cost = np.where(urban_mask, cost + 1000, cost) # zona urbana: casi prohibido
    return cost, s, d_infra

def route(cost, start_rc, end_rc):
    path, total = route_through_array(cost, start_rc, end_rc,
                                      fully_connected=True, geometric=True)
    return np.array(path), total

def straight_path(start_rc, end_rc):
    n = int(max(abs(end_rc[0]-start_rc[0]), abs(end_rc[1]-start_rc[1]))) + 1
    r = np.linspace(start_rc[0], end_rc[0], n).round().astype(int)
    c = np.linspace(start_rc[1], end_rc[1], n).round().astype(int)
    return np.stack([r, c], 1)

def metrics(path, dem, cost, d_infra, res, s, urban, water, major):
    step = np.hypot(*np.diff(path, axis=0).T) * res
    z = dem[path[:, 0], path[:, 1]]
    dz = np.diff(z)
    u = urban[path[1:, 0], path[1:, 1]]
    w = water[path[:, 0], path[:, 1]]
    cruces = int(w[0]) + int((w[1:] & ~w[:-1]).sum())
    mj = major[path[:, 0], path[:, 1]]
    cruces_ruta = int(mj[0]) + int((mj[1:] & ~mj[:-1]).sum())
    return {
        "longitud_km": round(step.sum() / 1000, 2),
        "desnivel_subida_m": round(float(dz[dz > 0].sum()), 1),
        "desnivel_bajada_m": round(float(-dz[dz < 0].sum()), 1),
        "pendiente_max_pct": round(float(s[path[:, 0], path[:, 1]].max()), 1),
        "pendiente_media_pct": round(float(s[path[:, 0], path[:, 1]].mean()), 1),
        "pct_paralelo_a_infraestructura": round(
            100 * float((d_infra[path[:, 0], path[:, 1]] <= NEAR_INFRA_M).mean()), 1),
        "m_en_zona_edificada": round(float(step[u].sum()), 0),
        "cruces_de_rutas_principales": cruces_ruta,
        "cruces_de_cursos_de_agua": cruces,
        "m_sobre_cursos_de_agua": round(float(step[w[1:]].sum()), 0),
        "costo_acumulado": round(float(cost[path[:, 0], path[:, 1]].sum()), 0),
    }

# ----------------------------- DATOS REALES -----------------------------
def load_real():
    import rasterio, geopandas as gpd, osmnx as ox
    from rasterio.windows import from_bounds
    from rasterio.features import rasterize
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_origin
    from pyproj import Transformer
    from shapely.geometry import Point

    tr = Transformer.from_crs("EPSG:4326", CRS_UTM, always_xy=True)
    ax, ay = tr.transform(A[1], A[0]); bx, by = tr.transform(B[1], B[0])
    xmin, xmax = min(ax, bx) - BUFFER_M, max(ax, bx) + BUFFER_M
    ymin, ymax = min(ay, by) - BUFFER_M, max(ay, by) + BUFFER_M
    xmin, ymin = np.floor(xmin/RES)*RES, np.floor(ymin/RES)*RES
    w, h = int((xmax-xmin)//RES)+1, int((ymax-ymin)//RES)+1
    dst_tf = from_origin(xmin, ymin + h*RES, RES, RES)

    # --- DEM Copernicus GLO-30 (tile 1°x1°; si el área cruza tiles, ver nota README)
    def tile_url(lat, lon):
        ns = "S" if lat < 0 else "N"; ew = "W" if lon < 0 else "E"
        la = int(np.floor(lat)); lo = int(np.floor(lon))
        name = f"Copernicus_DSM_COG_10_{ns}{abs(la):02d}_00_{ew}{abs(lo):03d}_00_DEM"
        return f"/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"
    os.environ["AWS_NO_SIGN_REQUEST"] = "YES"
    dem = np.zeros((h, w), dtype="float32")
    inv = Transformer.from_crs(CRS_UTM, "EPSG:4326", always_xy=True)
    corners = [inv.transform(x, y) for x in (xmin, xmax) for y in (ymin, ymax)]
    tiles = {(int(np.floor(la)), int(np.floor(lo))) for lo, la in corners}
    for la, lo in tiles:
        with rasterio.open(tile_url(la + .5, lo + .5)) as src:
            reproject(rasterio.band(src, 1), dem, dst_transform=dst_tf, dst_crs=CRS_UTM,
                      resampling=Resampling.bilinear, dst_nodata=np.nan)

    # --- OSM
    ll = inv.transform(xmin, ymin) + inv.transform(xmax, ymax)   # (lon0,lat0,lon1,lat1)
    bbox = (ll[0], ll[1], ll[2], ll[3])                           # left,bottom,right,top
    ENDPOINTS = ["https://overpass.kumi.systems/api",
                 "https://overpass.private.coffee/api",
                 "https://lz4.overpass-api.de/api",
                 "https://maps.mail.ru/osm/tools/overpass/api",
                 "https://overpass-api.de/api"] * 2
    import time
    def osm(tags, required=False):
        last = None
        for url in ENDPOINTS:
            for attr in ("overpass_url", "overpass_endpoint"):
                if hasattr(ox.settings, attr):
                    setattr(ox.settings, attr, url)
            ox.settings.requests_timeout = 60
            try:
                g = ox.features_from_bbox(bbox=bbox, tags=tags)
                print("OSM OK", url, list(tags)[0], len(g))
                return g.to_crs(CRS_UTM)
            except Exception as e:
                last = e
                msg = str(e)
                if "InsufficientResponse" in type(e).__name__ or "No matching features" in msg:
                    print("OSM sin elementos para", list(tags)[0]); 
                    return gpd.GeoDataFrame(geometry=[], crs=CRS_UTM)
                print("Falla", url, "->", msg[:90]); time.sleep(3)
        if required:
            raise RuntimeError("No se pudo descargar OSM (%s) de ningún servidor. "
                               "Reintentá en unos minutos: %s" % (list(tags)[0], last))
        return gpd.GeoDataFrame(geometry=[], crs=CRS_UTM)
    roads = osm({"highway": ["motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "track"]}, required=True)
    pipes = osm({"man_made": "pipeline"})
    water = osm({"waterway": True}, required=True)
    urban = osm({"landuse": ["residential", "commercial", "industrial"], "place": ["town", "village"]})
    bld = osm({"building": True})

    def rast(gdf, buf=0):
        if len(gdf) == 0:
            return np.zeros((h, w), bool)
        geoms = [g.buffer(buf) if buf else g for g in gdf.geometry if g is not None]
        return rasterize(geoms, out_shape=(h, w), transform=dst_tf, fill=0, default_value=1,
                         all_touched=True).astype(bool)
    PRINCIPALES = ["motorway", "trunk", "primary", "secondary", "tertiary"]
    major_gdf = roads[roads["highway"].isin(PRINCIPALES)] if len(roads) else roads
    infra = rast(roads, 5) | rast(pipes, 5)
    return dict(dem=dem, infra=infra, major=rast(major_gdf, 5), water=rast(water, 15), urban=rast(urban) | rast(bld, 50),
                tf=dst_tf, tr=tr, inv=inv, roads=roads, pipes=pipes, water_gdf=water,
                start=(ax, ay), end=(bx, by))

# ----------------------------- DATOS SINTÉTICOS (solo test) ---------------
def load_synthetic():
    rng = np.random.default_rng(0)
    h, w = 400, 600
    y, x = np.mgrid[0:h, 0:w]
    dem = 400 + 30*np.sin(x/90) + 20*np.cos(y/70) + 0.02*x
    dem += 80 * (np.abs(y - 200) < 10)            # barda (escarpa)
    dem[(np.abs(y-200) < 10) & (np.abs(x-420) < 25)] -= 80   # paso natural
    infra = np.zeros((h, w), bool); infra[290:293, :] = True
    water = np.zeros((h, w), bool); water[:, 300:303] = True
    urban = np.zeros((h, w), bool); urban[250:330, 100:200] = True
    return dict(dem=dem.astype("float32"), infra=infra, major=infra.copy(), water=water, urban=urban)

# ----------------------------- SALIDAS ----------------------------------
def save_png(dem, cost, p_opt, p_str, path_png):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    ls = LightSource(azdeg=315, altdeg=40)
    shade = ls.hillshade(dem, vert_exag=3, dx=RES, dy=RES)
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(shade, cmap="gray"); ax.imshow(dem, cmap="terrain", alpha=.45)
    ax.plot(p_str[:, 1], p_str[:, 0], "--", c="white", lw=1.5, label="Línea recta")
    ax.plot(p_opt[:, 1], p_opt[:, 0], c="red", lw=2.5, label="Traza optimizada")
    ax.scatter(*p_opt[0][::-1], c="blue", s=60, zorder=5, label="A")
    ax.scatter(*p_opt[-1][::-1], c="black", s=60, zorder=5, label="B")
    ax.legend(); ax.set_title("Traza optimizada vs línea recta")
    ax.axis("off"); fig.savefig(path_png, dpi=200, bbox_inches="tight"); plt.close(fig)

def save_anim(dem, infra, urban, p_opt, p_str, path_gif, path_mp4):
    """GIF (y MP4 si hay ffmpeg): primero la recta, después la traza optimizada."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter, writers
    from matplotlib.colors import LightSource, ListedColormap
    from matplotlib.patches import Patch
    h, w = dem.shape
    shade = LightSource(azdeg=315, altdeg=40).hillshade(dem, vert_exag=3, dx=RES, dy=RES)
    fig, ax = plt.subplots(figsize=(9, 9 * h / w))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.imshow(shade, cmap="gray"); ax.imshow(dem, cmap="terrain", alpha=.35)
    ax.imshow(np.ma.masked_where(~urban, urban.astype(float)), cmap=ListedColormap(["#ff8c00"]), alpha=.55)
    ax.imshow(np.ma.masked_where(~infra, infra.astype(float)), cmap=ListedColormap(["#222222"]), alpha=.8)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.axis("off")
    ax.scatter(p_opt[0][1], p_opt[0][0], c="blue", s=70, zorder=6)
    ax.scatter(p_opt[-1][1], p_opt[-1][0], c="black", s=70, zorder=6)
    ax.legend(handles=[Patch(color="#ff8c00", alpha=.6, label="Zona edificada"),
                       Patch(color="#222222", label="Caminos y ductos existentes")],
              loc="upper right", fontsize=8)
    l_str, = ax.plot([], [], "--", c="white", lw=2)
    l_opt, = ax.plot([], [], c="red", lw=3)
    txt = ax.text(0.02, 0.04, "", transform=ax.transAxes, fontsize=12, color="white",
                  bbox=dict(facecolor="black", alpha=.6, pad=5))
    def cum_km(p):
        return np.r_[0, np.cumsum(np.hypot(*np.diff(p, axis=0).T) * RES)] / 1000
    ks, ko = cum_km(p_str), cum_km(p_opt)
    F1, F2, HOLD = 45, 70, 30
    def update(f):
        if f < F1:
            k = max(2, int(len(p_str) * (f + 1) / F1))
            l_str.set_data(p_str[:k, 1], p_str[:k, 0])
            l_opt.set_data([], [])   # reinicia: el MP4 se guarda tras el GIF
            txt.set_text(f"Línea recta: {ks[k-1]:.1f} km")
        else:
            k = max(2, int(len(p_opt) * min(f - F1 + 1, F2) / F2))
            l_str.set_data(p_str[:, 1], p_str[:, 0])
            l_opt.set_data(p_opt[:k, 1], p_opt[:k, 0])
            txt.set_text(f"Traza optimizada: {ko[k-1]:.1f} km")
        return l_str, l_opt, txt
    anim = FuncAnimation(fig, update, frames=F1 + F2 + HOLD, blit=False)
    anim.save(path_gif, writer=PillowWriter(fps=20), dpi=100)
    if writers.is_available("ffmpeg"):
        anim.save(path_mp4, writer=FFMpegWriter(fps=20, bitrate=2500), dpi=100)
    plt.close(fig)

def save_folium(d, p_opt, p_str, path_html):
    import folium
    from rasterio.transform import xy
    def ll(p):
        xs, ys = xy(d["tf"], p[:, 0], p[:, 1])
        lon, lat = d["inv"].transform(np.array(xs), np.array(ys))
        return list(zip(lat, lon))
    m = folium.Map(location=[(A[0]+B[0])/2, (A[1]+B[1])/2], zoom_start=12, tiles="OpenStreetMap")
    folium.TileLayer("Esri.WorldImagery", name="Satélite").add_to(m)
    for g, col, nm in [(d["roads"], "gray", "Rutas/caminos"), (d["pipes"], "black", "Ductos OSM")]:
        if len(g): folium.GeoJson(g.to_crs(4326)[["geometry"]], name=nm,
                                   style_function=lambda f, c=col: {"color": c, "weight": 1.5}).add_to(m)
    folium.PolyLine(ll(p_str), color="blue", weight=2, dash_array="6", tooltip="Línea recta").add_to(m)
    folium.PolyLine(ll(p_opt), color="red", weight=4, tooltip="Traza optimizada").add_to(m)
    folium.Marker(A, tooltip="A: Loma Campana").add_to(m)
    folium.Marker(B, tooltip="B: Facilidad central (hipotética)").add_to(m)
    folium.LayerControl().add_to(m); m.save(path_html)

def main():
    os.makedirs(OUT, exist_ok=True)
    synthetic = "--synthetic" in sys.argv
    d = load_synthetic() if synthetic else load_real()
    dem = d["dem"]
    if synthetic:
        s_rc, e_rc = (40, 60), (360, 540)
    else:
        from rasterio.transform import rowcol
        s_rc = rowcol(d["tf"], *d["start"]); e_rc = rowcol(d["tf"], *d["end"])
    # Las instalaciones de origen y destino (edificios) no deben bloquear la entrada
    rr, cc = np.ogrid[:dem.shape[0], :dem.shape[1]]
    for (r0, c0) in (s_rc, e_rc):
        d["urban"][(rr - r0) ** 2 + (cc - c0) ** 2 <= (400 / RES) ** 2] = False
    cost, s, d_infra = build_cost(dem, RES, d["infra"], d["water"], d["urban"], d["major"])
    p_opt, _ = route(cost, s_rc, e_rc)
    p_str = straight_path(s_rc, e_rc)
    m_opt = metrics(p_opt, dem, cost, d_infra, RES, s, d["urban"], d["water"], d["major"])
    m_str = metrics(p_str, dem, cost, d_infra, RES, s, d["urban"], d["water"], d["major"])
    res = {"optimizada": m_opt, "linea_recta": m_str,
           "ahorro_costo_pct": round(100 * (1 - m_opt["costo_acumulado"] / m_str["costo_acumulado"]), 1)}
    json.dump(res, open(f"{OUT}/metricas.json", "w"), indent=2, ensure_ascii=False)
    print(json.dumps(res, indent=2, ensure_ascii=False))
    save_png(dem, cost, p_opt, p_str, f"{OUT}/traza.png")
    try:
        save_anim(dem, d["infra"], d["urban"], p_opt, p_str,
                  f"{OUT}/animacion_traza.gif", f"{OUT}/animacion_traza.mp4")
    except Exception as e:
        print("No se pudo generar la animación:", e)
    if not synthetic:
        save_folium(d, p_opt, p_str, f"{OUT}/mapa_interactivo.html")

if __name__ == "__main__":
    main()
