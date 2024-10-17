import pymatgen.core as mg
from multiprocessing import Pool
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import sys
# from IPython.display import Image, display
# %matplotlib inline

name = sys.argv[1]
structure = mg.Structure.from_str(open(name + ".cif").read(), fmt="cif")
a, b, c = structure.lattice.abc
alpha, beta, gamma, = structure.lattice.angles
print (a, b, c)
print (alpha, beta, gamma)
elements = [x.specie.name for x in structure.sites]
coords = [x.frac_coords for x in structure.sites]

xrd_c = XRDCalculator( wavelength="CuKa1")

def compare_XRD(a_scale, b_scale, c_scale, alpha_scale, beta_scale, gamma_scale):
    new_lattice = mg.Lattice.from_parameters(a * a_scale, b* b_scale, c* c_scale, alpha *alpha_scale , beta * beta_scale, gamma * gamma_scale)
    new = mg.Structure(new_lattice, elements, coords)
#     new.to(filename="new.cif")
   
    # c.show_plot(structure)
    xrd = xrd_c.get_pattern(new, two_theta_range=(0, 26.5))

    target_100 = 5.1
    target_010 = 5.1
    target_200 = 10.2
    target_020 = 10.2
    target_001 = 25.7
    tolerance = 0.5
    condition010 = False
    condition100 = False
    condition200 = False
    condition020 = False
    condition001 = False
    inten200 = 0
    inten020 = 0
    inten100 = 0
    inten010 = 0
#     condition = False
    for two_theta, i, hkls, d_hkl in zip(xrd.x, xrd.y, xrd.hkls, xrd.d_hkls):
        if hkls[0]["hkl"] == (1, 0, 0):
            cal_100 = two_theta
            inten100 = i
            if abs(target_100 - cal_100)<=tolerance:
                condition100 = True
            else:
                condition100 = False

        if hkls[0]["hkl"] == (0, 1, 0):
            cal_010 =  two_theta
            inten010 = i
            if abs(target_010 - cal_010)<=tolerance:
                condition010 = True
            else:
                condition010 = False

        if hkls[0]["hkl"] == (2, 0, 0):
            cal_200 = two_theta
            inten200 = i
            if abs(target_200 - cal_200)<=tolerance:
                condition200 = True
            else:
                condition200 = False

        if hkls[0]["hkl"] == (0, 2, 0):
            cal_020 = two_theta
            inten020 = i
            if abs(target_020 - cal_020)<=tolerance:
                condition020 = True
            else:
                condition020 = False

        if hkls[0]["hkl"] == (0, 0, 1):
            cal_001 =  two_theta
            if abs(target_001 - cal_001)<=tolerance:
                condition001 = True
            else:
                condition001 = False

    if inten100 > inten010:
        if condition100 == True:
            condition010 = True
        elif condition100 == False:
            condition010 = False
    if inten100 < inten010:
        if condition010 == True:
            condition100 = True
        elif condition010 == False:
            condition100 = False  

    if inten200 > inten020:
        if condition200 == True:
            condition020 = True
        elif condition200 == False:
            condition020 = False
    if inten200 < inten020:
        if condition020 == True:
            condition200 = True
        elif condition020 == False:
            condition200 = False    

    condition = condition100 and condition010 and condition020 and condition200 and condition001
 
    if condition == True:
    	print (a_scale, b_scale, c_scale, alpha_scale, beta_scale, gamma_scale,condition)
    	filename = "cif/"+str(a_scale)+"_"+ str(b_scale) +"_"+ str(c_scale) +"_"+ str(alpha_scale) +"_"+ str(beta_scale) +"_"+ str(gamma_scale)
    	new.to(filename=filename + ".cif")


    return  [a_scale, b_scale, c_scale, alpha_scale, beta_scale, gamma_scale, condition]

if __name__ == "__main__":
	pool = Pool(processes=1000)
	for a_scale in range(90,110,2):
	    for b_scale in range(90,110,2):
	        for c_scale in range(90,110,2):
	            for alpha_scale in range(90,110,2):
	                for beta_scale in range(90,110,2):
	                    for gamma_scale in range(90,110,2):
	                        pool.apply(compare_XRD, args = (a_scale/100, b_scale/100, c_scale/100, alpha_scale/100, beta_scale/100, gamma_scale/100, ))
	                        #data = compare_XRD(a_scale/100, b_scale/100, c_scale/100, alpha_scale/100, beta_scale/100, gamma_scale/100)
                        # if data[-1] == True:
                        # 	print (data)
	print("start")
	pool.close()
	pool.join()



# a_scale = 100
# b_scale = 100
# c_scale = 100
# alpha_scale = 100
# beta_scale  = 100
# gamma_scale = 100
# data = compare_XRD(a_scale/100, b_scale/100, c_scale/100, alpha_scale/100, beta_scale/100, gamma_scale/100)
# print(data)
