import trio
import asks
import time

async def geobarometer(session, id):
    url = "http://geobarometer.default.enki-portal.net"
    params = {
        'SiO2': "77.5", 
        'TiO2':  "0.08", 
        'Al2O3':"12.5", 
        'Fe2O3': "0.207",
        'Cr2O3': "0.0", 
        'FeO':   "0.473", 
        'MnO':   "0.0",
        'MgO':   "0.03", 
        'NiO':   "0.0", 
        'CoO':   "0.0",
        'CaO':   "0.43", 
        'Na2O':  "3.98", 
        'K2O':   "4.88", 
        'P2O5':  "0.0", 
        'H2O':  "10.0",
        'fO2 offset':"0.0"
    }
    print (id, 'launched.')
    response = await session.get(url, params=params)
    print (id, response.content)
    return response.content

async def run_sequence(number):
    from asks.sessions import Session
    session = Session("http://geobarometer.default.enki-portal.net", connections=1)
    async with trio.open_nursery() as n:
        for id in range(0,number):
            n.start_soon(geobarometer, session, id)
ncase = 100

time_ini = time.perf_counter()
trio.run(run_sequence, ncase)
time_fin = time.perf_counter()

print ("Time elapsed:", time_fin - time_ini, "Average:", (time_fin - time_ini)/ncase)

