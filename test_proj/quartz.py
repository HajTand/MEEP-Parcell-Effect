import math

import matplotlib.pyplot as plt
import numpy as np
from meep.materials import fused_quartz
"из библиотеки материалов берем квартц (со своим eps = eps (w) экспериментальным)"

import meep as mp

"количество микселей на микрометров"
"размер ячейки 5 нм = 1/200"
resolution = 200  # pixels/μm

"dpml - граничное условие, чем больше, тем лучше поглощает (не отражает)"
"d - значит слой в микрометрах"
dpml = 1.0
"sz - размер рабочей области целого вместе с pml"
sz = 10 + 2 * dpml
cell_size = mp.Vector3(z=sz)
"размер ячейки"
pml_layers = [mp.PML(dpml)]
"количество слоев pml???"

wvl_min = 0.4
wvl_max = 0.8
fmin = 1 / wvl_max
fmax = 1 / wvl_min
fcen = 0.5 * (fmax + fmin)
df = fmax - fmin
nfreq = 50

"задаем источник и его параметры, центр геометрический"
"ширины df источник"
sources = [ 
    mp.Source(
        mp.GaussianSource(fcen, fwidth=df),
        component=mp.Ex,
        center=mp.Vector3(z=-0.5 * sz + dpml),
    )
]
"задаем симуляцию, dim - размерность моделирования"
"размер ячейки - cell_size, исчтоник - sources, разрешение - resolution"
"boundary layers - ???"
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    sources=sources,
    dimensions=1,
    resolution=resolution,
)
"задаем монитор, поток"
refl_fr = mp.FluxRegion(center=mp.Vector3(z=-0.25 * sz))
refl = sim.add_flux(fcen, df, nfreq, refl_fr)
"окно df количеством nfreq"
"монитор где записываются данные"

"запуск симуляции, задаем когда закончится симуляция"
"каждые 50 фемтосекунд после того как закончился импульс света"
"компонента Ex должно стать 1e-9"
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-9))

"создание экрана"
empty_flux = mp.get_fluxes(refl)
empty_data = sim.get_flux_data(refl)
sim.reset_meep()
"запуск заного с материалами"

"добавляем тела внуть рабочей области"
"материал наш - fused_quartz"
geometry = [
    mp.Block(
        mp.Vector3(mp.inf, mp.inf, 0.5 * sz),
        center=mp.Vector3(z=0.25 * sz),
        material=fused_quartz,
    )
]

sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    dimensions=1,
    resolution=resolution,
)

refl = sim.add_flux(fcen, df, nfreq, refl_fr)
sim.load_minus_flux_data(refl, empty_data)

"sim.run - запуск симуляции"
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(), 1e-9))

"поток flux = get_fluxes - считывание с монитора"
refl_flux = mp.get_fluxes(refl)
R_meep = -1 * np.divide(refl_flux, empty_flux)

freqs = mp.get_flux_freqs(refl)
wvls = np.divide(1, freqs)

eps_quartz = (
    lambda l: 1
    + 0.6961663 * math.pow(l, 2) / (pow(l, 2) - pow(0.0684043, 2))
    + 0.4079426 * pow(l, 2) / (pow(l, 2) - pow(0.1162414, 2))
    + 0.8974794 * pow(l, 2) / (pow(l, 2) - pow(9.896161, 2))
)
R_fresnel = lambda l: math.pow(
    math.fabs(1 - math.sqrt(eps_quartz(l))) / (1 + math.sqrt(eps_quartz(l))), 2
)
R_analytic = [R_fresnel(i) for i in wvls]

plt.figure()
plt.plot(wvls, R_meep, "bo-", label="meep")
plt.plot(wvls, R_analytic, "rs-", label="analytic")
plt.xlabel("wavelength (μm)")
plt.ylabel("reflectance")
plt.axis([0.4, 0.8, 0.0340, 0.0365])
plt.xticks(list(np.arange(0.4, 0.9, 0.1)))
plt.legend(loc="upper right")
plt.show()
