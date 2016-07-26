import csv
import sys

try: paraview.simple
except: from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

f = open('DeformedGeo_RBFref.csv', 'rt')
index = 0
points = []
try:
    reader = csv.reader(f)
    for row in reader:
	if index == 0:
            index = index + 1
        elif index > 0:
            Sphere1 = Sphere()
            Sphere1.Radius = 0.05
            Sphere1.Center = [float(row[0]), float(row[1]), float(row[2])]
            Sphere1.ThetaResolution = 10
            Sphere1.PhiResolution = 10
            RenderView = GetRenderView()
            DataRepresentation = Show()
            DataRepresentation.DiffuseColor = [1.0, 0.0, 0.0]
finally:
    f.close()

Render()


