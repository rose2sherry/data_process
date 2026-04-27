import ee
import geopandas as gpd
import json
import os

# ==============================================
# 0. 授权与初始化
# ==============================================
PROJECT_ID = 'ee-studyuse188124'
try:
    # 尝试直接使用项目 ID 初始化
    ee.Initialize(project=PROJECT_ID)
except Exception as e:
    # 如果没授权过，先授权
    print("正在尝试重新授权...")
    ee.Authenticate()
    ee.Initialize(project=PROJECT_ID)

print(f"✅ GEE 初始化成功！当前项目: {PROJECT_ID}")
print("开始处理任务...")

# ==============================================
# 1. 本地路径与配置区
# ==============================================
# 本地 SHP 文件的完整路径
LOCAL_SHP_PATH = r"G:\全国实验\356个有城市边界数据\杭州市.shp"

# 导出到 Google Drive 的配置
DRIVE_FOLDER = 'gee_python'
FILE_NAME = 'hangzhou'

# ==============================================
# 2. 读取本地 SHP 并强制转码及提取边界框
# ==============================================
if not os.path.exists(LOCAL_SHP_PATH):
    print(f"❌ 找不到文件，请确认路径是否正确: {LOCAL_SHP_PATH}")
    exit()

try:
    print(f"⏳ 正在读取本地 SHP: {os.path.basename(LOCAL_SHP_PATH)}")
    gdf = gpd.read_file(LOCAL_SHP_PATH)

    if gdf.crs is None:
        print("⚠️ 警告: 该 SHP 没有定义坐标系，假设为 WGS84。")
    elif gdf.crs.to_string() != "EPSG:4326":
        print(f"🔄 检测到坐标系 {gdf.crs.to_string()}，正在转为 EPSG:4326 (WGS84)...")
        gdf = gdf.to_crs("EPSG:4326")

    # 获取该 SHP 的绝对经纬度极值 (minx, miny, maxx, maxy)
    minx, miny, maxx, maxy = gdf.total_bounds
    print(f"📍 经纬度极值: 经度({minx:.4f} ~ {maxx:.4f}), 纬度({miny:.4f} ~ {maxy:.4f})")

    # 依然保留原始的多边形几何体，用于最后 GEE 的精细筛选
    ee_features = []
    for _, row in gdf.iterrows():
        geom_dict = row['geometry'].__geo_interface__
        ee_features.append(ee.Feature(ee.Geometry(geom_dict)))

    roi_geometry = ee.FeatureCollection(ee_features).geometry()
    print("✅ 边界几何数据准备就绪。")

except Exception as e:
    print(f"❌ 处理本地文件出错: {e}")
    exit()

# ==============================================
# 3. 本地 Python 预计算相交瓦片 (物理隔绝法)
# ==============================================
tile_paths = [
    # 北纬 50° ~ 55°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n55_e120_n50",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n55_e125_n50",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e125_n55_e130_n50",
    # 北纬 45° ~ 50°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e080_n50_e085_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e085_n50_e090_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e090_n50_e095_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n50_e115_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n50_e120_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n50_e125_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e125_n50_e130_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e130_n50_e135_n45",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e135_n50_e140_n45",
    # 北纬 40° ~ 45°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e070_n45_e075_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e075_n45_e080_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e080_n45_e085_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e085_n45_e090_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e090_n45_e095_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e095_n45_e100_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e100_n45_e105_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n45_e110_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n45_e115_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n45_e120_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n45_e125_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e125_n45_e130_n40",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e130_n45_e135_n40",
    # 北纬 35° ~ 40°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e070_n40_e075_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e075_n40_e080_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e080_n40_e085_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e085_n40_e090_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e090_n40_e095_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e095_n40_e100_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e100_n40_e105_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n40_e110_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n40_e115_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n40_e120_n35",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n40_e125_n35",
    # 北纬 30° ~ 35°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e075_n35_e080_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e080_n35_e085_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e085_n35_e090_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e090_n35_e095_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e095_n35_e100_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e100_n35_e105_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n35_e110_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n35_e115_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n35_e120_n30",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n35_e125_n30",
    # 北纬 25° ~ 30°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e085_n30_e090_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e090_n30_e095_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e095_n30_e100_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e100_n30_e105_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n30_e110_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n30_e115_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n30_e120_n25",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n30_e125_n25",
    # 北纬 20° ~ 25°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e095_n25_e100_n20",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e100_n25_e105_n20",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n25_e110_n20",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n25_e115_n20",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n25_e120_n20",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e120_n25_e125_n20",
    # 北纬 15° ~ 20°
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e105_n20_e110_n15",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e110_n20_e115_n15",
    "projects/sat-io/open-datasets/GLOBAL_BUILDING_ATLAS/e115_n20_e120_n15"
]

print("⏳ 正在本地计算需要加载的瓦片，彻底抛弃无关区域...")
valid_paths = []

for path in tile_paths:
    # 提取最后的名字，例如 e120_n35_e125_n30
    basename = path.split('/')[-1]
    parts = basename.split('_')

    if len(parts) == 4:
        # 解析该瓦片自带的经纬度范围 (忽略字母 e 和 n)
        lon1 = float(parts[0][1:])
        lat1 = float(parts[1][1:])
        lon2 = float(parts[2][1:])
        lat2 = float(parts[3][1:])

        t_minx = min(lon1, lon2)
        t_maxx = max(lon1, lon2)
        t_miny = min(lat1, lat2)
        t_maxy = max(lat1, lat2)

        # 判断瓦片边界是否与我们提取的 minx, miny, maxx, maxy 有相交
        if (minx <= t_maxx and maxx >= t_minx and miny <= t_maxy and maxy >= t_miny):
            valid_paths.append(path)

if not valid_paths:
    print("❌ 警告：您的本地 SHP 范围没有覆盖到任何内置的建筑瓦片！")
    exit()

print(f"🎯 Python 预处理完成！成功锁定 {len(valid_paths)} 个真实相交瓦片:")
for vp in valid_paths:
    print(f"   👉 {vp.split('/')[-1]}")

# ==============================================
# 4. GEE 精确合并与导出
# ==============================================
print("⏳ 正在 GEE 中合并目标瓦片并进行精确筛选...")

# 因为我们已经极大地缩减了瓦片数量（通常只剩 1 到 2 个），
# 这里改用更稳定的 .merge() 来代替容易引起索引崩溃的 .flatten()
clipped_buildings = ee.FeatureCollection(valid_paths[0])
for i in range(1, len(valid_paths)):
    clipped_buildings = clipped_buildings.merge(ee.FeatureCollection(valid_paths[i]))

# 使用真实的城市多边形边界进行最后一次精确过滤
clipped_buildings = clipped_buildings.filterBounds(roi_geometry)

# 附加清洗：确保要素有高度属性，防止空边界混入
clipped_buildings = clipped_buildings.filter(ee.Filter.notNull(['height']))

print("🚀 正在提交导出任务到 Google Drive...")
task = ee.batch.Export.table.toDrive(
    collection=clipped_buildings,
    description=f'Task_{FILE_NAME}',
    folder=DRIVE_FOLDER,
    fileNamePrefix=FILE_NAME,
    fileFormat='SHP'
)

task.start()

print("-" * 30)
print(f"✨ 任务已成功提交！")
print(f"📂 云盘文件夹: {DRIVE_FOLDER}")
print(f"📄 文件名: {FILE_NAME}.zip")
print("-" * 30)