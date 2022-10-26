#!/usr/bin/env python
import string,sys,os,copy

"""
Evaluate operator strings, and generate tex file
"""

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, Dagger)
from sympy import (
    Add, Mul, symbols, expand, pprint, Rational, latex, Dummy, KroneckerDelta
)


#latex equation output for "short" expressions
def Equation2Tex(Equ,file):

    file.write(r'\begin{eqnarray} \begin{aligned}')
    file.write("\n")
    file.write(Equ+r'\\')
    file.write("\n")
    file.write("\end{aligned} \end{eqnarray} \n")

#latex equation output for "long" expressions
def Equation2Tex_Resize(Equ,file):

    file.write(r'\begin{eqnarray} \begin{aligned}')
    file.write("\n")
    file.write(r'\resizebox{.8\hsize}{!}{$'+Equ+r'$}\\')
    file.write("\n")
    file.write("\end{aligned} \end{eqnarray} \n")

# For stringsssss!
def Permute_Op(Opstring):

    result    = []
    for str in Opstring.args:
        print(str)
       
        c_part     = []
        string1    = []
        for factor in str.args:
            if factor.is_commutative:
                c_part.append(factor)
            else:
                string1.append(factor)

        #initialize the sign
        stringlist = []
        sign       = []
        sign.append(1)
        stringlist.append(string1)
        i=0
        normal_ordering(stringlist, sign, i)
        i=0
        for factor in stringlist:
            result.append(Mul(*c_part)*sign[i]*Mul(*factor))                    
            i=i+1
    print("final results:")
    print(result)
    return Add(*result)

#For string!
def Permute_Str(Opstring):
    result     = []   
    c_part     = []
    string1    = []
    for factor in Opstring.args:
        if factor.is_commutative:
            c_part.append(factor)
        else:
            string1.append(factor)
    print(c_part)
    print(string1)

    #initialize the sign
    stringlist = []
    sign       = []
    sign.append(1)
    stringlist.append(string1)
    i=0
    normal_ordering(stringlist, sign,i)
    i=0
    for factor in stringlist:
        result.append(Mul(*c_part)*sign[i]*Mul(*factor))                    
        i=i+1
    print("final results:")
    print(result)
    return Add(*result)

# Judge the status of a string
def is_normal_ordered(string1):

    NOP=[]
    NCreator=0
    NAnnihal=0
    n = len(string1)

    for i in range(len(string1)):
        if string1[i].is_commutative:
            pass
        else:
            if string1[i].is_only_q_creator:
                NOP.append(1)
                NCreator=NCreator+1
            elif string1[i].is_only_q_annihilator:
                NOP.append(-1)
                NAnnihal=NAnnihal+1
            else:
                NOP.append(0)

    #check creators
    N=len(NOP)
    NC=0
    for i in range(0,N):
       if NOP[i]==1:
           NC=NC+1
       elif NOP[i]!=1:
           break
    #print "testing NC %d" %NC
    NA=0   
    for j in range(N-1,-1,-1):
        if NOP[j]==-1:
            NA=NA+1
        elif NOP[j]!=-1:
            break
    #print "testing NA %d" %NA    
    if (NA==NAnnihal) & (NC==NCreator):
        print(string1)
        print ("This string contains %d creators and %d annihilators, which is normal ordered" %(NCreator, NAnnihal))
        return True
    else:
        print(string1)
        print ("This string contains %d creators and %d annihilators, which is NOT! normal ordered" %(NCreator, NAnnihal))
        return False

#Judge the right hand side of a string, no "creators" on the right side.
def locus_of_right_creator(string1):

    NOP=[]
    NCreator=0
    n = len(string1)
    for i in range(0,n):
        if string1[i].is_commutative:
            pass
        else:
            if string1[i].is_only_q_creator:
                NOP.append(i)
                NCreator=NCreator+1
            else:
                NOP.append(-1)

    #check creators
    NC=0
    N=0
    for i in range(n):
       if NOP[i]!=-1:
           NC=NC+1
       else:
           break
    if NC==NCreator:
        #print(string1)
        #print "All %d creators is normal ordered" %(NCreator)
        return -1
    else:
        #print(string1)
        j=0
        for i in range(n): 
            if NOP[i]!=-1:
                j=j+1
                if j==NC+1:
                    break
            else:
                pass
        N=NOP[i]
        #print "The creators on %d is not normal ordered" %(N)
        return N

#Judge the right hand side of a string, no "annihilators" on the left side.
def locus_of_left_annihilator(string1):

    NOP=[]
    NAnnihilator=0
    n = len(string1)
    for i in range(0,n):
        if string1[i].is_commutative:
            pass
        else:
            if string1[i].is_only_q_annihilator:
                NOP.append(i)
                NAnnihilator=NAnnihilator+1
            else:
                NOP.append(-1)

    #check annihilators
    NA=0
    m=len(NOP)
    for i in range(m-1,-1,-1):
       if NOP[i]!=-1:
           NA=NA+1
       else:
           break

    if NA==NAnnihilator:
        #print "All %d annihilators is normal ordered" %(NAnnihilator)
        return -1
    else:
        j=0
        for i in range(m-1,-1,-1):
            if NOP[i]!=-1:
                j=j+1
                if j==NA+1:
                    break
            else:
                pass
        N=NOP[i]
        #print "The annihilator on %d is not normal ordered" %(N)
        return N


#Normal ordering strings recursively!!
#'i' is the layer of strings,
# it will run over 'all' strings, including newly generated strings.
#'results' contains normal ordered strings
#'sign' contains sign for each term in 'results'
def normal_ordering(results, sign, i):

    # resutls is a list of string
    NS = len(results)
    #print i, NS
    if i == NS:
        return
    #print results[i]
    if is_normal_ordered(results[i]):
        i=i+1
        normal_ordering(results, sign, i)
        return 
    # Move creator on the right to left 
    if locus_of_right_creator(results[i])==-1:
        # Move annihilater on the left to the right
        if locus_of_left_annihilator(results[i])==-1:
            print ("Error !!!!!!")
        else:
            n = locus_of_left_annihilator(results[i])
            Delta1 = (results[i][n+1].op_symbol=='f') & (results[i][n].op_symbol=='f+')
            Delta2 = (results[i][n+1].op_symbol=='f+') & (results[i][n].op_symbol=='f')
            print (Delta1,Delta2, n)
            if Delta1 | Delta2:
                # add new string to both resutls and sign
                tmp = results[i][:]
                tmp.append(KroneckerDelta(tmp[n+1].state,tmp[n].state))
                tmp.pop(n+1)
                tmp.pop(n)
                if tmp:
                    results.append(tmp)
                    tmp = sign[i]
                    sign.append(tmp)
                # switch elements
                results[i][n+1], results[i][n] = results[i][n], results[i][n+1]
                sign[i]=sign[i]*-1
            else:
                # switch elements
                results[i][n+1], results[i][n] = results[i][n], results[i][n+1]
                sign[i]=sign[i]*-1
    else:
        n = locus_of_right_creator(results[i]) 
        Delta1 = (results[i][n-1].op_symbol=='f') & (results[i][n].op_symbol=='f+')
        Delta2 = (results[i][n-1].op_symbol=='f+') & (results[i][n].op_symbol=='f')
        if Delta1 | Delta2:
            # add new string to both resutls and sign
            tmp = results[i][:]
            tmp.append(KroneckerDelta(tmp[n-1].state,tmp[n].state))
            tmp.pop(n)
            tmp.pop(n-1)
            if tmp:
                results.append(tmp)
                tmp = sign[i]
                sign.append(tmp)
            # switch elements
            results[i][n-1], results[i][n] = results[i][n], results[i][n-1]
            sign[i]=sign[i]*-1
            #print(sign[i],results[i])
            #normal_ordering(results, sign, i)
        else:
            # switch elements
            results[i][n-1], results[i][n] = results[i][n], results[i][n-1]
            sign[i]=sign[i]*-1
            #print(sign[i],results[i])
            #normal_ordering(results, sign, i)
     
    normal_ordering(results, sign, i)    

# Generate RPA equations

def main():

    base     = sys.argv[0]
    filename = os.path.splitext(base)[0]
    fo = open(filename+'.tex', "w")

    fo.write(r'\documentclass{article}')
    fo.write("\n")
    fo.write(r'\usepackage[fleqn]{amsmath}')
    fo.write("\n")
    fo.write(r'\usepackage{latexsym, amsfonts}')
    fo.write("\n")
    fo.write(r'\usepackage{pdflscape}')
    fo.write("\n")
    fo.write(r'\parindent=0pt\relax')
    fo.write("\n")
    fo.write(r'\begin{document}')
    fo.write("\n")
    fo.write(r'\begin{landscape}')
    fo.write("\n")

    # Hamiltonian
    p, q, r, s = symbols('p,q,r,s', cls=Dummy)
    f = AntiSymmetricTensor('f', (p,), (q,))
    pr = Fd(p)*F(q)
    v = AntiSymmetricTensor('v', (p, q), (r, s))
    pqsr = Fd(p)*Fd(q)*F(s)*F(r)

    H = f*pr + Rational(1, 4)*v*pqsr
    # Print Hamiltonian
    fo.write("The spin orbital hamiltonian:\n")
    Equation2Tex(latex(H),fo)

    # Exciation and de-excitation operators
    i, j = symbols('i,j', below_fermi=True)
    a, b = symbols('a,b', above_fermi=True)
    ia = Fd(a)*F(i)
    jb = Fd(j)*F(b)

    C = Commutator
    #First com
    comm1= C(H,ia)
    comm2= C(jb,comm1)
    fo.write("Evaluating operators:\n")
    Equation2Tex(latex(comm2),fo)

    comm2= comm2.doit()
    comm2= comm2.expand()
    #comm2= Permute_Str(comm2)
    comm2= Permute_Op(comm2)

    Equation2Tex(latex(comm2),fo) 

    #finish writing
    fo.write("\end{landscape}\n")
    fo.write("\end{document}\n")
    fo.close()
    os.system("pdflatex "+filename+'.tex')

if __name__ == "__main__":
    main()
