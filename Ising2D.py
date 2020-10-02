import numpy
from collections import defaultdict
from matplotlib import pyplot
import itertools
#%matplotlib inline

J = 1*(1.38e-23)
kB = 1.38e-23
red = 10  

sitios = list()
spins = dict()
vecinos = defaultdict(list)

for x, y in itertools.product(range(red), range(red)):
    sitios.append((x,y)) 

print(sitios)

def configuración_aleatoria_spin():
    for spin in sitios:
        spins[spin] = numpy.random.choice([-1, 1])
configuración_aleatoria_spin()
#print(sites)


vecinos = defaultdict(list)
for site in spins:
    x, y = site
    if x + 1 < red:
        vecinos[site].append(((x + 1) % red, y))
    if x - 1 >= 0:
        vecinos[site].append(((x - 1) % red, y))
    if y + 1 < red:
        vecinos[site].append((x, (y + 1) % red))
    if y - 1 >= 0:
        vecinos[site].append((x, (y - 1) % red))

def energía_sitio(site):
    Energía = 0.0
    for nbh in vecinos[site]:
        Energía += spins[site] * spins[nbh]
    return -J * Energía

def energía_total():
    Energía = 0.0
    for site in sitios:
        Energía += energía_sitio(site)
    return 0.5 * Energía

def magnetización():
    mag = 0.0
    for spin in spins.values():
        mag += spin
    return mag

#plot_spins()
print("Magnetización = ", magnetización())
#magnetization =  0.0

def metropolis(site, T):
    Spin_f = spins[site]
    Energía_final = energía_sitio(site)
    spins[site] *= -1
    Energía_inicial = energía_sitio(site)
    deltaE = Energía_inicial - Energía_final
    if deltaE <= 0:
        pass
    else:
        if numpy.random.uniform(0, 1) <= numpy.exp(-deltaE/(kB*T)):
            pass
        else:
            spins[site] *= -1

def monte_carlo(T):
    for i in range(len(sitios)):
        int_rand_site = numpy.random.randint(0, len(sitios))
        rand_site = sitios[int_rand_site]
        metropolis(rand_site, T)

pasos_simulación = 4000
T_f= 15.0
T_i = 0.01
stepT = 0.1


temps = numpy.arange(T_i, T_f, stepT)
energías = numpy.zeros(shape=(len(temps), pasos_simulación))
magnetizaciones = numpy.zeros(shape=(len(temps), pasos_simulación))
configuración_aleatoria_spin()
for ind_T, T in enumerate(temps):
    for i in range(pasos_simulación):
        monte_carlo(T)
        energías[ind_T, i] = energía_total()
        magnetizaciones[ind_T, i] = magnetización()

tau = pasos_simulación // 2
energía_media = numpy.mean(energías[:, tau:], axis=1)
magnetización_media = abs(numpy.mean(magnetizaciones[:, tau:], axis=1))

print(energía_media)
print("Energía media = ",sum(energía_media)//pasos_simulación)
print(magnetización_media)

pyplot.figure()
pyplot.plot(energía_media, label="Energía", linestyle=':', color='DarkBlue', marker='o')
pyplot.legend()
pyplot.xlabel(r"tiempo ($s$)")
pyplot.ylabel(r"Energía ($J$)")
pyplot.grid()
pyplot.show()

pyplot.figure()
pyplot.plot(temps, magnetización_media, label="Magnetización", linestyle=':', color='DarkBlue', marker='o')
pyplot.legend()
pyplot.xlabel(r"Temperatura ($K$)")
pyplot.ylabel(r"Magnetización ($A/m$)")
pyplot.grid()
pyplot.show()




