import random
import time


def is_prime(n, k=5):
    """Miller-Rabin algorithm"""
    if n <= 1:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # Write n as 2^r * d + 1
    r = 0
    d = n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a


def multiplicative_inverse(e, phi):
    """Calculate the multiplicative inverse of e mod phi."""
    d = 0
    x1, x2 = 0, 1
    y1, y2 = 1, 0
    temp_phi = phi

    while e > 0:
        temp1 = temp_phi // e
        temp2 = temp_phi - temp1 * e
        temp_phi = e
        e = temp2

        x = x2 - temp1 * x1
        y = y2 - temp1 * y1

        x2 = x1
        x1 = x
        y2 = y1
        y1 = y

    if temp_phi == 1:
        d = y2 + phi

    return d


def generate_pq(keysize):
    prime = random.randrange(2 ** (keysize - 1), 2 ** keysize)
    while not is_prime(prime):
        prime = random.randrange(2 ** (keysize - 1), 2 ** keysize)
    return prime


def generate_keypair(keysize):

    p=0
    q=0
    while not is_prime(p):
        p = random.randrange(2 ** (keysize - 1), 2 ** keysize)

    while not is_prime(q):
        q = random.randrange(2 ** (keysize - 1), 2 ** keysize)

    n = p * q
    phi = (p - 1) * (q - 1)

    e = random.randrange(1, phi)
    g = gcd(e, phi)
    while g != 1:
        e = random.randrange(1, phi)
        g = gcd(e, phi)

    d = multiplicative_inverse(e, phi)

    return (e, n), (d, n)


def encrypt_rsa(message, public_key):
    e, n = public_key
    cipher = [pow(ord(char), e, n) for char in message]
    return cipher


def decrypt_rsa(cipher, private_key):
    d, n = private_key
    message = [chr(pow(char, d, n)) for char in cipher]
    return ''.join(message)


message = "UZAYCETINKAYA"

print("Generating RSA key ...")
public_key, private_key = generate_keypair(1024)
print("\nPublic Key (e, n):", public_key)
print("Private Key (d, n):", private_key)
print("\nEncrypting message...")
start_time = time.time()
encrypted_message = encrypt_rsa(message, public_key)
end_time = time.time()
print("Encryption time:", end_time - start_time, "seconds")
print("Encrypted Message:", encrypted_message)

print("\nDecrypting message...")
start_time = time.time()
decrypted_message = decrypt_rsa(encrypted_message, private_key)
end_time = time.time()
print("Decryption time:", end_time - start_time, "seconds")
print("Decrypted Message:", decrypted_message)


# S-DES ALGORITHM
# --------- KEY GENERATION PART ---------

def sdes_randomkey():
    key = ""
    for _ in range(10):
        key += str(random.randint(0, 1))
    return key


def p10_function(key):
    # P10 permutation table from S-DES.pdf
    P10 = [3, 5, 2, 7, 4, 10, 1, 9, 8, 6]

    permuted_key = [key[i - 1] for i in P10]
    return permuted_key

def p8_function(key):
    # P8 permutation table
    P8 = [6, 3, 7, 4, 8, 5, 10, 9]
    # Permute the key using P8 table
    permuted_key = [key[i - 1] for i in P8]
    return permuted_key

def circular_left_shift(bits, shift):

    return bits[shift:] + bits[:shift]


key = "1010000010"  # 10-bit key from sdes_randomkey function
print("Original Key:", key)

# Permute key
p10_key = p10_function(key)


left_half = p10_key[:5]
right_half = p10_key[5:]


left_half_shifted = circular_left_shift(left_half, 1)
right_half_shifted = circular_left_shift(right_half, 1)
shifted_key = left_half_shifted + right_half_shifted

k1 = p8_function(shifted_key)

left_half_shifted_2 = circular_left_shift(left_half_shifted, 2)
right_half_shifted_2 = circular_left_shift(right_half_shifted, 2)
shifted_key_2 = left_half_shifted_2 + right_half_shifted_2

k2 = p8_function(shifted_key_2)


print("K1:", k1)
print("K2:", k2)

# --------- ENCRYPTION PART ---------




#def encrypt_sdes():


#def decrypt_sdes():

'''
# encryption

def des_encrypt(plaintext, key):
    state = plaintext

    # subkey generation (scheduling) = generate 2 subkeys
    subkeys = generate_subkey(key)

    # initial permutation
    state = ip(state)

    # 2 fiestel networks
    for i in range(2):
        state = feistel(state, subkeys[i])

    # final permutation
    state = ip_1(state)
    
    ciphertext = state

    return ciphertext


def generate_subkey(key):
    # P10
    key = p10(key)

    # split key
    key_l = key[:5]
    key_r = key[5:]
    
    # LS-1
    key_l_shifted = left_shift(key_l, 1)
    key_r_shifted = left_shift(key_r, 1)

    key = key_l_shifted + key_r_shifted

    # P8 --> subkey 0
    subkey_0 = p8(key)

    # LS-2 
    key_l_shifted = left_shift(key_l_shifted, 2)
    key_r_shifted = left_shift(key_r_shifted, 2)

    key = key_l_shifted + key_r_shifted

    # P8 --> subkey 1
    subkey_1 = p8(key)

    return subkey_0, subkey_1


def left_shift(input, shift_amount):
    return input[shift_amount:] + input[:shift_amount]


def feistel(state, subkey):
    # split state
    state_l = state[:4]
    state_r = state[4:]

    # first round function
    f_state_r = round_function(state_r, subkey)

    # XOR right side with left side of the state
    state_l = xor(state_l, f_state_r)

    state = state_r + state_l

    return state


def round_function(state, subkey):
    # E/P (Expansion Permutation)
    state = e_p(state)
    
    # XOR with subkey
    state = xor(state, subkey)

    # Sbox
    state = sbox(state)

    # P4
    state = p4(state)

    return state


# === permutation & substitution functions === #
def p10(input):
    P10 = [2, 4, 1, 6, 3, 9, 0, 8, 7, 5]  # Fixed P10 permutation
    output = [input[i] for i in P10]
    return ''.join(output)
    
def p8(input):
    P8 = [5, 2, 6, 3, 7, 4, 9, 8]  # Fixed P8 permutation
    output = [input[i] for i in P8]
    return ''.join(output)

def p4(input):
    P4 = [1, 3, 2, 0]  # Fixed P4 permutation
    output = [input[i] for i in P4]
    return ''.join(output)

# IP (Initial Permuation)
def ip(input):
    IP = [1, 5, 2, 0, 3, 7, 4, 6]  # Fixed IP permutation
    output = [input[i] for i in IP]
    return ''.join(output)

# IP**-1 (Final Permutation)
def ip_1(input):
    IP_1 = [3, 0, 2, 4, 6, 1, 7, 5]  # Fixed IP-1 permutation
    output = [input[i] for i in IP_1]
    return ''.join(output)

# E/P (Expansion Permutation)
def e_p(input):
    E_P = [3, 0, 1, 2, 1, 2, 3, 0]  # Fixed E/P permutation
    output = [input[i] for i in E_P]
    return ''.join(output)

def sbox(input):
    S0 = [['01', '00', '11', '10'],
          ['11', '10', '01', '00'],
          ['00', '10', '01', '11'],
          ['11', '01', '11', '10']]
    
    S1 = [['00', '01', '10', '11'],
          ['10', '00', '01', '11'],
          ['11', '00', '01', '00'],
          ['10', '01', '00', '11']]

    input_0 = input[:4] 
    input_1 = input[4:]

    row_s0 = int(input_0[0] + input_0[3], 2)
    col_s0 = int(input_0[1] + input_0[2], 2)
    row_s1 = int(input_1[0] + input_1[3], 2)
    col_s1 = int(input_1[1] + input_1[2], 2)

    output = S0[row_s0][col_s0] + S1[row_s1][col_s1]
    return output

def xor(a, b):
    result = ''
    for i in range(len(a)):
        result += str(int(a[i]) ^ int(b[i]))
    return result

print("---------------------------S-DES Algorithm---------------------------")

key = '1010000010'
plaintext = '11010111'

print("Key:", key)
print("Plaintext:", plaintext)

ciphertext = des_encrypt(plaintext, key)
print("Ciphertext:", ciphertext)

'''


