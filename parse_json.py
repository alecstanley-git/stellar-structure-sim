import json

with open("out.json", "r") as f:
    data = json.load(f)

snap0 = data[0]
snap10 = data[-1]

print("AGE 0:")
print(f"Pc={snap0['pc_cgs']:.3e} Tc={snap0['tc_k']:.3e} R={snap0['r_rsun']:.3f}")
print("AGE 10:")
print(f"Pc={snap10['pc_cgs']:.3e} Tc={snap10['tc_k']:.3e} R={snap10['r_rsun']:.3f}")
