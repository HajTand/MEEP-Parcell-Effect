susceptibilities = [mp.LorentzianSusceptibility(frequency=1.1, gamma=1e-5, sigma=0.5),
                    mp.LorentzianSusceptibility(frequency=0.5, gamma=0.1, sigma=2e-5)]
"восприимчивость хи - поляризуемость на единицу объема"
"поляризуемость это про одну молекулу"


default_material = mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)
"epsilon 2.25 - постоянный вклад"
"гамма - частота соударения ( вещество с поглощением )"
"freq = w0"