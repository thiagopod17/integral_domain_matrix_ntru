######################################
# Algebra Coefficients Z(sqrt(-3))
######################################
# define Ring
K.<t> = NumberField(x^2+3)

# modular operation as defined in the paper
def z_mod_p(z,p):
    a = int(z[0])
    b = int(z[1])
    return K([a%p,b%p])


def center_lift_coef(a,q):
    a_real = a[0]
    a_img = a[1]
    if a_real > q/2:
        a_real -= q
    if a_img > q/2:
        a_img -= q
    return K([a_real,a_img])

######################################
# Matrices Coefficients Z(sqrt(-3))
######################################
# Add, Multiply works using A+B and A*B directly

# modular matrix operations
def matrixMODp(A,p):
    return matrix([[z_mod_p(ele,p) for ele in l] for l in A])


# generate a random element 
import numpy as np
def generate_random_element(N, d):
    """
    Generate a random element from {-1, 0, 1} with the following probabilities:
      - P(-1) = d / (2*N**2)
      - P( 1) = d / (2*N**2)
      - P( 0) = (2*N**2 - 2*d) / (2*N**2)
    """
    # Calculate probabilities
    p_minus1 = d / (2*N**2)
    p_zero   = (2*N**2 - 2*d) / (2*N**2)
    p_plus1  = d / (2*N**2)
    
    # Generate and return the random element
    return int(np.random.choice([-1, 0, 1], p=[p_minus1, p_zero, p_plus1]))

# generate a random matrix with elements in K
def random_matrix(N, d):
  # Check if the norm is 0 modulo p, in which case the inverse does not exist.
  if 2*N*N - 2*d < 0:
    raise ValueError(" 2*N*N - 2*d < 0")
  return matrix([[K([generate_random_element(N, d),generate_random_element(N, d)]) for i in range(N)] for j in range(N)])


# generate random matrix with fixed numbers of +1 and -1 
def random_matrix_fixed_d(N, d):
    for trial in range(1000):
        matrix_rnd = random_matrix(N, d)
        count_one = 0
        count_minus_one = 0
        for line in matrix_rnd:
            for element in line:
                for val in range(2):
                    if element[val] == 1:
                        count_one += 1
                    if element[val] == -1:
                        count_minus_one += 1
        if count_one == d and count_minus_one == d:
            return matrix_rnd


def random_matrix_with_integer_determinant_coprime_p_q(N,p,q,d):
    """
    Return a matrix of dimension NxN whose determinant is non zero, integer and 
    it is coprime with both p and q, which are coprime integers as well. 
    This returned matrix can be used as the private key of the matrix ntru algorithm 

    Input:
    N,p,q,d: matrix ntru parameters

    Output:
    found: boolean indicating whether a random matrix was found satisfying the desired properties
    F: if found is true, the matrix F. Otherwise, return false. 
    """
    found = False
    for i in range(1000):
        F = random_matrix_fixed_d(N,d)
        F_det = F.determinant()
        if F_det in ZZ and F_det != 0:
            F_det_int = int(F_det[0])
            if abs(gcd(F_det_int,p)) == 1:
                if abs(gcd(F_det_int,q)) == 1:
                    found = True
                    return found,F
    return found,found



def key_gen_F_with_integer_determinant(N,p,q,d):
    """
    Key generation function for the integral domain matrix ntru, where F has integer 
    determinant that is coprime with p and q. 

    Input:
    N,p,q,d: matrix ntru parameters

    Output:
    F,Fp,Fq,G,H: matrices representing the system private and public keys
    """
    _,F = random_matrix_with_integer_determinant_coprime_p_q(N,p,q,d)
    det_F = F.determinant()
    gcd_p,alpha_p,beta_p = xgcd(det_F,p)
    gcd_q,alpha_q,beta_q = xgcd(det_F,q)
    F_adjugate = F.adjugate()
    Fp = matrixMODp(alpha_p * F_adjugate,p)
    Fq = matrixMODp(alpha_q * F_adjugate,q)
    G = random_matrix_fixed_d(N,d)
    H = matrixMODp(p * Fq * G, q)
    return F,Fp,Fq,G,H











# ----------------------------------------------------------------
# Matrices for N = 3, crafted to satisfy properties and to go to the boundaries of decryption failure to test
# robustness of the theorem 1 or spot errors of the original paper
def random_matrix_F_with_N_3_d_5_fixed_line():
    found = False
    for i in range(1000):
        F_start_matrix = random_matrix_with_fixed_line()
        F_det = F_start_matrix.determinant()
        if F_det in ZZ and F_det != 0:
            F_det_int = int(F_det[0])
            if abs(gcd(F_det_int,p)) == 1:
                if abs(gcd(F_det_int,q)) == 1:
                    found = True
                    return found, matrix(F_start_matrix)


def random_matrix_with_fixed_line():
    A = [[K([-1,-1]), K([1,-1]), K([-1,-1])],
    [K([0,0]), K([0,0]), K([0,0])],
    [K([0,0]), K([0,0]), K([0,0])]]
    post_to_choose = [ [1,0,0], [1,0,1], [1,1,0], [1,1,1], [1,2,0], [1,2,1], [2,0,0], [2,0,1], [2,1,0], [2,1,1], [2,2,0], [2,2,1] ]
    index_choosen = np.random.choice(range(len(post_to_choose)), size = 4, replace = False)  
    for pos in index_choosen:
        matrix_pos = post_to_choose[pos]
        if matrix_pos[2]:
            A[matrix_pos[0]][matrix_pos[1]] += K([1,0])
        else:
            A[matrix_pos[0]][matrix_pos[1]] += K([0,1])
    return matrix(A)


def random_matrix_with_fixed_column():
    A = [[K([1,-1]), K([0,0]), K([0,0])],
    [K([-1,-1]), K([0,0]), K([0,0])],
    [K([1,-1]), K([0,0]), K([0,0])]]
    post_to_choose = [ [0,1,0], [0,1,1], [0,2,0], [0,2,1], [1,1,0], [1,1,1], [1,2,0], [1,2,1], [2,1,0], [2,1,1], [2,2,0], [2,2,1] ]
    index_choosen = np.random.choice(range(len(post_to_choose)), size = 4, replace = False)
    cont = 0  
    for pos in index_choosen:

        # first coef is the -1 position, the others are +1 coefs.
        if cont == 0:
            change_sign = -1
        else:
            change_sign = 1
        cont += 1

        # get matrix position and put a +1 or -1 
        matrix_pos = post_to_choose[pos]
        if matrix_pos[2] == 1:
            A[matrix_pos[0]][matrix_pos[1]] += change_sign * K([1,0])
        else:
            A[matrix_pos[0]][matrix_pos[1]] += change_sign * K([0,1])
    return matrix(A)





# generate random matrix with fixed numbers of +1 and -1 
def count_ones_and_minus_ones(A):
    count_one = 0
    count_minus_one = 0
    for line in A:
        for element in line:
            for val in range(2):
                if element[val] == 1:
                    count_one += 1
                if element[val] == -1:
                    count_minus_one += 1
    return count_one, count_minus_one

    
def random_matrix_fixed_d_fast(N,d):
    if d > N*N:
        raise ValueError("d > N*N")
    A = [[K([0,0])]*N for i in range(N)]
    a = np.random.choice(range(2*N*N),size = 2*d,replace=False)



    a_positive = a[0:d]
    a_negative = a[d:2*d]


    for i in a_positive:
        if i < N*N:
            idiv = i // N
            imod = i % N
            A[idiv][imod] += K([1,0])
        else:
            idiv = (i-N*N) // N
            imod = (i-N*N) % N
            A[idiv][imod] += K([0,1])


    for i in a_negative:
        if i < N*N:
            idiv = i // N
            imod = i % N
            A[idiv][imod] += K([-1,0])
        else:
            idiv = (i-N*N) // N
            imod = (i-N*N) % N
            A[idiv][imod] += K([0,-1])
    return matrix(A)





######################################
# Encrypt and Decrypt
######################################
def encrypt(M,H,N,p,q):
    """
    Encrypt message M given public key H

    Input:
    N,p,q: dimension of matrix and primes p and q
    M,H: message and public key H
    Output:
    E: encrypted message
    """
    R = random_matrix_fixed_d_fast(N, d)
    E = H * R + M

    return matrixMODp(E,q)

def decrypt(E,F,F_inv_p,N,p,q):
    """
    Decrypt encrypted message E given private key F and its inverse mod p

    Input:
    E: encrypted message
    F,F_inv_p: public key F and its inverse mod p

    Output:
    M: message
    """
    A = matrixMODp(F * E,q)
    # Choose coefficients of A into −q/2 to q/2.
    A_center_lift = matrix([[center_lift_coef(ele,q) for ele in l] for l in A])
    B = matrixMODp(A_center_lift,p)
    C = matrixMODp(F_inv_p * B,p)
    C_center_lift = matrix([[center_lift_coef(ele,p) for ele in l] for l in C])
    return C_center_lift

######################################
# New key gen, no restriction on integer determinant or alpha or beta integer in the diofantine equation
######################################
#-------------------------
# modular reduction
#-------------------------
def aMODb(a, b, returnQuotientAndRest = 0):
    N = b.norm(); u = a*b.conjugate()

    qa,ra = divmod(ZZ(u[0]),ZZ(N)) # the argument here is u =  a*b.conjugate()
    qb,rb = divmod(ZZ(u[1]),ZZ(N))
    if(ra > N/2):
        qa = qa + 1
    if(rb > N/2):
        qb = qb + 1
    q =  K([qa,qb])
    r =  a - q*b
    if(returnQuotientAndRest):
        return q,r
    else: 
        return r


#-------------------------
# extendedEuclidean
#-------------------------
def extendedEuclidean(f,g):
    # changed the notation to f and g 
    # to make it easier to follow the algorithm on p. 48 of Garthen and Gerhard book.
    # returns
    # s,t,r such that s * f + t * g == r
    r0 = f
    s0 = K([1,0])
    t0 = K([0,0])
    r1 = g
    s1 = K([0,0])
    t1 = K([1,0])
    i = 1
    ri = r1
    rnext = r1
    rprev = r0
    sprev = s0
    si = s1
    tprev = t0
    ti = t1
    while (ri != K([0,0]) ):
        qi,rnext = aMODb(rprev,ri,1)
        snext = sprev - qi * si
        tnext = tprev - qi * ti
        i = i + 1
        
        tprev = ti
        ti = tnext
        
        sprev = si
        si = snext  
        
        rprev = ri
        ri = rnext
    return sprev,tprev,rprev


def keygen_new(N,p,q,d,maxit = 100):
    key_found = False
    for i in range(maxit):
        A = random_matrix_fixed_d_fast(N,d)
        A_det = A.determinant()
        alpha_p,beta_p,gcd_det_p = extendedEuclidean(A_det,K([p,0]))
        if gcd_det_p == 1 or gcd_det_p == -1:
            alpha_q,beta_q,gcd_det_q = extendedEuclidean(A_det,K([q,0]))
            if gcd_det_q == 1 or gcd_det_q == -1:
                A_inv_p = matrixMODp(gcd_det_p * alpha_p * A.adjugate(),p)
                A_inv_q = matrixMODp(gcd_det_q * alpha_q * A.adjugate(),q)
                key_found == True
                G = random_matrix(N, d)
                H = matrixMODp(p * A_inv_q * G,q)
                return A,A_inv_p,A_inv_q,G,H
    return 'failed'




# decryption failure condition get smallest q value 
# q has a minimum value so that there is zero decryption failure 
# this function get a prime q satisfying it
prime_list = list(primes(10000))
def q_min_value(N,p,prime_list = prime_list):
    lim_inf = 8*N*(p + ceil((p-1)/2) )   
    for number in prime_list:
        if number > lim_inf:
            val = number
            break
    return val