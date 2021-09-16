#Run with sage
import os

#Find out the name of the directory containing this file
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

R.<a,b,c,d,e,f,g,h> = PolynomialRing(QQ, 8);

def poly_C_code(poly_str, #string to be converted
                var_list, #list of variables in poly_str
                acbs, #corresponding names of acb elements
                filename, #output file name (in this directory)
                res = "res", #name of output variable
                term = "temp", #name of variable storing terms
                prec = "prec", #name of prec variable
                pow_suffix = "_pows" #use acb_ptrs named var1_pows, etc
):
    n = len(var_list)
    assert (n <= 8) #if not, use another poly ring...

    #Convert poly_str to a,b,c,d,e,f
    #Here we need the names of the variables to not overlap.
    for k in [0..n-1]:
        var = var_list[k]
        new_var = ["a","b","c","d","e","f","g","h"][k]
        poly_str = poly_str.replace(var, new_var)
    #print("Poly: ", poly_str)

    #Build polynomial
    poly = R(poly_str);
    #Compute partial degrees
    degrees = []
    for k in [0..n-1]:
        deg = poly.degree([a,b,c,d,e,f,g,h][k])
        degrees += [deg]
    #print("Partial degrees:", degrees)

    #Make list of power list names
    power_lists = [var_list[k] + pow_suffix for k in [0..n-1]]
    
    complete_filename = os.path.join(__location__, filename)

    #It is up to the user to ensure that all acb variables are
    #correctly declared; we simply init/clear them here
    with open(complete_filename, 'a') as file:
        file.write("\n")
        
        file.write("/* Init all power lists */\n")
        for k in [k for k in [0..n-1] if (degrees[k] > 0)]:
            file.write(power_lists[k] + " = _acb_vec_init({});\n".format(degrees[k]+1))
        file.write("\n")

        file.write("/* Precompute all power lists */\n")
        for k in [k for k in [0..n-1] if (degrees[k] > 0)]:
            file.write("acb_one(&{}[0]);\n".format(power_lists[k]))
            for l in [1..degrees[k]]:
                file.write("acb_mul(&{vec}[{nextindex}], &{vec}[{oldindex}], {acb}, {prec});\n".format(vec=power_lists[k], nextindex=l, oldindex=l-1, acb=acbs[k], prec=prec))
            file.write("\n")

        file.write("/* Add successive terms to res */\n")
        file.write("acb_zero({});\n".format(res))

        for m in [0..len(poly.monomials())-1]:
            monomial = poly.monomials()[m]
            coef = poly.coefficients()[m]
            exponents = monomial.exponents()[0] #This is a tuple of exponents
            
            file.write("\n")
            file.write("acb_one({});\n".format(term))
            for k in [0..n-1]:
                if exponents[k] != 0:
                    file.write("acb_mul({term}, {term}, &{powlist}[{exp}], {prec});\n".format(term=term, powlist=power_lists[k], exp=exponents[k], prec=prec))
            file.write("acb_mul_si({term}, {term}, {coef}, {prec});\n".format(term=term, coef=coef.numerator(), prec=prec))
            file.write("acb_div_si({term}, {term}, {coef}, {prec});\n".format(term=term, coef=coef.denominator(), prec=prec))
            file.write("acb_add({res}, {res}, {term}, {prec});\n".format(res=res, term=term, prec=prec))
        file.write("\n")

        file.write("/* Clear all power lists */\n")
        for k in [k for k in [0..n-1]if (degrees[k] > 0)]:
            file.write("_acb_vec_clear({}, {});\n".format(power_lists[k], degrees[k]+1))
        
        file.write("\n")

def declare_power_lists(acbs, # names of acb elements
                        filename,
                        pow_suffix = "_pows"):
    n = len(acbs)
    assert (n <= 8) #if not, use another poly ring...
                     
    #Make list of power list names
    power_lists = [acbs[k] + pow_suffix for k in [0..n-1]]
    
    complete_filename = os.path.join(__location__, filename)

    with open(complete_filename, 'a') as file:
        for k in [0..n-1]:
            file.write("acb_ptr {};\n".format(power_lists[k]))

def make_helper_function(inputfile, #string to be converted
                var_list, #list of variables in poly_str
                acbs, #corresponding names of acb elements
                function_name, #output file name (in this directory)
                res = "res", #name of output variable
                term = "temp", #name of variable storing terms
                prec = "prec", #name of prec variable
                pow_suffix = "_pows" #use acb_ptrs named var1_pows, etc
):
    with open(os.path.join(__location__, inputfile), 'r') as file:
        poly_str = file.read().replace('\n', '')
    filename = function_name + ".c"
    complete_filename = os.path.join(__location__, filename)
    n = len(var_list)

    #Define function
    with open(complete_filename, 'w') as file: #Erase content
        file.write("\nvoid " + function_name + "(acb_t " + res + ", ")
        for k in [0..n-1]:
            file.write("const acb_t " + acbs[k] + ", ")
        file.write("slong " + prec + ")\n{\n")
        file.write("acb_t " + term + ";\n")
    #Declare power lists
    declare_power_lists(acbs, filename, pow_suffix)
    #Init temp
    with open(complete_filename, 'a') as file:
        file.write("acb_init({});\n".format(term))
    #Do the computations
    poly_C_code(poly_str, var_list, acbs, filename, res, term, prec)
    #Clear temp and close
    with open(complete_filename, 'a') as file:
        file.write("acb_clear({});\n".format(term))
        file.write("}\n")
    #Indent correctly (taken from Flint code conventions)
    os.system("indent -bap -blf -bli0 -cbi0 -cdw -cli4 -cs -i4 -l79 -nbad -nbc -nce -npcs -nprs -nut -pmt -psl -saf -sai -saw -sbi0 -ss -ts4 " + complete_filename)

a_list = ["a0", "a1", "a2", "a3", "a4", "a5", "a6"]
make_helper_function("I2", a_list, a_list, "igusa_I2_autogen")
make_helper_function("I4", a_list, a_list, "igusa_I4_autogen")
make_helper_function("I6prime", a_list, a_list, "igusa_I6prime_autogen")
make_helper_function("I10", a_list, a_list, "igusa_I10_autogen")

    
# clebsch_from_igusa
varlist_clebsch = ["A", "B", "C", "I2", "I4", "I6", "I10"]
polA_clebsch = "-I2/120"
polB_clebsch = "(I4+720*A^2)/6750"
polC_clebsch = "(I6-8640*A^3+108000*A*B)/202500"
polD_clebsch = "(I10+62208*A^5-972000*A^3*B-1620000*A^2*C+3037500*A*B^2+6075000*B*C)/(-4556250)"
filename_clebsch = "clebsch_from_igusa.c"
acbs_clebsch = varlist_clebsch

with open(os.path.join(__location__, filename_clebsch), 'w') as file:
    pass #erase
declare_power_lists(acbs_clebsch, filename_clebsch);
poly_C_code(polA_clebsch, varlist_clebsch, acbs_clebsch, filename_clebsch, "A")
poly_C_code(polB_clebsch, varlist_clebsch, acbs_clebsch, filename_clebsch, "B")
poly_C_code(polC_clebsch, varlist_clebsch, acbs_clebsch, filename_clebsch, "C")
poly_C_code(polD_clebsch, varlist_clebsch, acbs_clebsch, filename_clebsch, "D")

# R^2
varlist_R = ["I2", "I4", "I6", "I10"]
make_helper_function("R2", varlist_R, varlist_R, "igusa_R2_autogen")

# mestre_conic
varlist_conic_1 = ["A", "B", "C", "D"];
varlist_conic_2 = ["A11", "A22", "A33", "A23", "A31", "A12", "U", "DP"]
acbs_conic_1 = ["&ABCD[0]", "&ABCD[1]", "&ABCD[2]", "&ABCD[3]"]
acbs_conic_2 = ["&A_coefs[0]", "&A_coefs[1]", "&A_coefs[2]", "&A_coefs[3]",
          "&A_coefs[4]", "&A_coefs[5]", "U", "I10"];
filename_conic = "mestre_conic.c"

with open(os.path.join(__location__, filename_conic), 'w') as file:
    pass #erase

declare_power_lists(varlist_conic_1, filename_conic);
declare_power_lists(varlist_conic_2, filename_conic);

polys_conic_1 = ["2*C+A*B/3", "D", "B*D/2+2*C*(B^2+A*C)/9",
                 "B*(B^2+A*C)/3+C*(2*C+A*B/3)/3", "D", "2*(B^2+A*C)/3"]
for k in [0..5]:
    poly_C_code(polys_conic_1[k], varlist_conic_1, acbs_conic_1, filename_conic, acbs_conic_2[k])

polys_conic_2 = [
    "A11*U^2*DP^8", "A22*DP^10", "A33*U^8",
    "A23*DP^5*U^4", "A31*DP^4*U^5", "A12*U*DP^9"]

# Do this by hand.
# for k in [0..5]:
#     poly_C_code(polys_conic_2[k], varlist_conic_2, acbs_conic_2,
#                 filename_conic, "&C_coefs[{}]".format(k))
          

# mestre_cubic: similarly, do second step by hand
varlist_cubic = ["A", "B", "C", "D"];
acbs_cubic_1 = ["&ABCD[0]", "&ABCD[1]", "&ABCD[2]", "&ABCD[3]"]
acbs_cubic_2 = ["&A_coefs[{}]".format(k) for k in [0..9]]
filename_cubic = "mestre_cubic.c"

with open(os.path.join(__location__, filename_cubic), 'w') as file:
    pass #erase

declare_power_lists(varlist_conic_1, filename_cubic);

polys_cubic = [
    "8*(A^2*C-6*B*C+9*D)/36",
    "4*(2*B^3+4*A*B*C+12*C^2+3*A*D)/36",
    "4*(A*B^3+4*A^2*B*C/3+4*B^2*C+6*A*C^2+3*B*D)/36",
    "4*(A*B^3+4*A^2*B*C/3+4*B^2*C+6*A*C^2+3*B*D)/36",
    "2*(2*B^4+4*A*B^2*C+4*A^2*C^2/3+4*B*C^2+3*A*B*D+12*C*D)/36",
    "2*(A*B^4+4*A^2*B^2*C/3+16*B^3*C/3+26*A*B*C^2/3+ 8*C^3+3*B^2*D+2*A*C*D)/36",
    "4*(3*B^4+6*A*B^2*C+8*A^2*C^2/3+2*B*C^2-3*C*D)/36",
    "2*(-2*B^3*C/3-4*A*B*C^2/3-4*C^3+9*B^2*D+8*A*C*D)/36",
    "2*(B^5+2*A*B^3*C+8*A^2*B*C^2/9+2*B^2*C^2/3 -B*C*D+9*D^2)/36",
    "1*(-2*B^4*C-4*A*B^2*C^2-16*A^2*C^3/9-4*B*C^3/3 +9*B^3*D+12*A*B*C*D+20*C^2*D)/36"]
for k in [0..9]:
    poly_C_code(polys_cubic[k], varlist_cubic, acbs_cubic_1, filename_cubic, acbs_cubic_2[k])
