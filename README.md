#  matrix NTRU over the integral domain 

### Implementation in sage of the system presented in the original paper and assessement of:
- Key generation
- Correct encryption and decryption
- Decryption failure conditions
- Generating counterexamples for original theorem
- Testing the corrected version of the theorem (to be published in SBSEG 2025)

### Low level functions include:
- module p operations over matrices 
- generation of random matrices satysfing specific formats
- extended euclidean algorithms following standard literature (see [2]).


### References
[1] 2023, IEEE, ICoCICs, Indah E. Wijayanti, Uha Isnaini, Anny Kartika Sari. Matrix NTRU Cryptosystem over Integral Domain from Published at IEEE - ICoCICs in 2023.
[2] von zur Gathen, J., & Gerhard, J. (2013). Modern Computer Algebra (3rd ed.). Cambridge University Press.

