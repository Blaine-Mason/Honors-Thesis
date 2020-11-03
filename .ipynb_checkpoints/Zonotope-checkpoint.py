from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import geopandas as gpd

p1 = Polygon([(0,0),(1,0),(1,1),(0,1)])
p2 = Polygon([(0,1),(1,1),(2,1),(2,2)])
newp = p1.union(p2)

p = gpd.GeoSeries(newp)
p.plot()
plt.show()