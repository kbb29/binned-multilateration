import chevron, pathlib, os, random

template = (pathlib.Path(__file__).parent / 'map.mustache').open('r').read()


def generateDonut(lat, lon, outerradius, innerradius, options={}) :
    return f'''L.donut([{lat}, {lon}], {{
        color: '{options.get('color', 'red')}',
        fillColor: '{options.get('color', 'red')}',
        fillOpacity: 0.4,
        radius: {outerradius},
        innerRadius: {innerradius}
    }})'''

def generateCircle(lat, lon, radius, options={}):
    return f'''L.circle([{lat}, {lon}], {{
        color: '{options.get('color', "red")}',
        fillColor: '{options.get('color', "red")}',
        fillOpacity: {options.get('opacity', "0.0")},
        radius: {radius}
    }})'''

def generateMarker(lat, lon, title):
    return f"L.marker([{lat}, {lon}], {{title: '{title}'}})"

def generateAddToFeaturesAndMap(featureCode):
    return f"features.push({featureCode}.addTo(map));\n"

def plotBeacons(beacons, preds=[], centroid=None, actual=None, options={}):
    featureJS = ""
    for i,beacon in enumerate(beacons):
        featureJS += generateAddToFeaturesAndMap(generateDonut(beacon.point.lat, beacon.point.lon, beacon.limits[1], beacon.limits[0], {'color': "red"}))
        featureJS += generateAddToFeaturesAndMap(generateMarker(beacon.point.lat, beacon.point.lon, f"kevin{i}"))

    if actual:
        featureJS += generateAddToFeaturesAndMap(generateMarker(actual.lat, actual.lon, 'actual'))
    
    for i,pred in enumerate(preds):
        featureJS += generateAddToFeaturesAndMap(generateMarker(pred.lat, pred.lon, f'pred{i}'))

    if centroid:
        featureJS += generateAddToFeaturesAndMap(generateCircle(centroid.point.lat, centroid.point.lon, centroid.radius, {'color': 'blue', 'opacity': "0.1"}))

    output = chevron.render(template, {'features': featureJS});
    
    tag = options.get('tag', str(random.random()))

    fn = pathlib.Path('/tmp') / f'map-{tag}.html'
    print('fn', fn)
    with fn.open('w') as fh:
        fh.write(output)

    print('fn2', fn)
    os.system(f'xdg-open {fn}')
