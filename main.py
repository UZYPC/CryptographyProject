import random
import time


# Function to check if a number is prime using Miller-Rabin algorithm
def is_prime(n, k=5):
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


# Function to calculate the greatest common divisor (GCD) of two numbers
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a


# Function to calculate the multiplicative inverse using Extended Euclidean algorithm
def multiplicative_inverse(e, phi):
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


# Function to generate a prime number of given key size
def generate_pq(keysize):
    prime = random.randrange(2 ** (keysize - 1), 2 ** keysize)
    while not is_prime(prime):
        prime = random.randrange(2 ** (keysize - 1), 2 ** keysize)
    return prime


# Function to generate public and private keys for RSA encryption
def generate_keypair(keysize):
    p = 0
    q = 0
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


# Function to encrypt a message using RSA algorithm
def encrypt_rsa(message, public_key):
    e, n = public_key
    cipher = [pow(ord(char), e, n) for char in message]
    return cipher


# Function to decrypt a message using RSA algorithm
def decrypt_rsa(cipher, private_key):
    d, n = private_key
    message = [chr(pow(char, d, n)) for char in cipher]
    return ''.join(message)


# Sample 10 messages (inputs) for RSA encryption
message1 = "UZAYCETINKAYA"
message2 = "OZANBASOGUL"
message3 = "DENIZOZKAHRAMAN"
message4 = "HELLO"
message5 = "WORLD"
message6 = "RSA"
message7 = "ENCRYPTION"
message8 = "DECRYPTION"
message9 = "ALGORITHM"
message10 = "SECURITY"

# Perform RSA encryption and decryption
print("Generating RSA key ...")

public_key, private_key = generate_keypair(1024)
print("\nPublic Key (e, n):", public_key)
print("Private Key (d, n):", private_key)

for i in range(1, 11):
    print("\nMessage", i, ":", eval("message" + str(i)))

    # Measure the time taken for encryption
    start_time = time.time()
    encrypted_message = encrypt_rsa(eval("message" + str(i)), public_key)
    end_time = time.time()
    print("Encryption time:", end_time - start_time, "seconds")
    print("Encrypted Message:", encrypted_message)

    print("\nMessage", i, "is now on decryption process...")
    # Measure the time taken for decryption
    start_time = time.time()
    decrypted_message = decrypt_rsa(encrypted_message, private_key)
    end_time = time.time()
    print("Decryption time:", end_time - start_time, "seconds")
    print("Decrypted Message:", decrypted_message)

print("\n")
print(
    "-----------------------------------------------------------------------------------------------------------------")
print("\n")

# S-DES ALGORITHM

# 10 plaintexts for testing and key for S-DES
plaintext1 = "11111111"
plaintext2 = "10101010"
plaintext3 = "11001100"
plaintext4 = "11110000"
plaintext5 = "00001111"
plaintext6 = "10110101"
plaintext7 = "01010101"
plaintext8 = "11100000"
plaintext9 = "00000000"
plaintext10 = "10010101"

# Shared (Secret) 10-bit key for S-DES
key = "1010000010"
print("Original Key:", key)

# S-DES S-Box matrices
S0 = [
    ["01", "00", "11", "10"],
    ["11", "10", "01", "00"],
    ["00", "10", "01", "11"],
    ["11", "01", "11", "10"]
]

S1 = [
    ["00", "01", "10", "11"],
    ["10", "00", "01", "11"],
    ["11", "00", "01", "00"],
    ["10", "01", "00", "11"]
]


# Function to convert binary string to decimal
def binary_to_decimal(binary_str):
    return int(binary_str, 2)


# Function to merge indices for S-Box lookup
def merge_indices(array):
    first_pair = ''.join([str(array[0]), str(array[3])])
    second_pair = ''.join([str(array[1]), str(array[2])])
    return first_pair, second_pair


# Function to perform S-Box lookup
def s_box_lookup(s_box, row, column):
    return s_box[row][column]


# Function to switch two bits (SW function in S-DES)
def switch(bit1, bit2):
    return bit2, bit1


# Function to generate a random 10-bit key for S-DES
def sdes_randomkey():
    key = ""
    for _ in range(10):
        key += str(random.randint(0, 1))
    return key


# Permutation functions for S-DES keys
def p10_function(key):
    P10 = [3, 5, 2, 7, 4, 10, 1, 9, 8, 6]

    permuted_key = [key[i - 1] for i in P10]
    return permuted_key


def ip_function(key):
    IP = [2, 6, 3, 1, 4, 8, 5, 7]

    permuted_key = [key[i - 1] for i in IP]
    return permuted_key


def ip_inv_function(key):
    IPinv = [4, 1, 3, 5, 7, 2, 8, 6]

    permuted_key = [key[i - 1] for i in IPinv]
    return permuted_key


def p8_function(key):
    P8 = [6, 3, 7, 4, 8, 5, 10, 9]

    permuted_key = [key[i - 1] for i in P8]
    return permuted_key


def p4_function(key):
    P4 = [2, 4, 3, 1]

    permuted_key = [key[i - 1] for i in P4]
    return permuted_key


# Key generation and permutation for S-DES
def circular_left_shift(bits, shift):
    return bits[shift:] + bits[:shift]


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


# Function to perform F-function and generate subkey
def F_R_subkey(R, subkey):
    subkey_first_half = subkey[:4]
    subkey_second_half = subkey[4:]

    def ep1_function(R):
        EP_1 = [4, 1, 2, 3]

        permuted_key = [R[i - 1] for i in EP_1]
        return permuted_key

    def ep2_function(R):
        EP_2 = [2, 3, 4, 1]

        permuted_key = [R[i - 1] for i in EP_2]
        return permuted_key

    p_array_1 = (ep1_function(R))
    p_array_2 = (ep2_function(R))

    first_half = [int(subkey_bit) ^ int(R_bit) for R_bit, subkey_bit in zip(p_array_1, subkey_first_half)]

    second_half = [int(subkey_bit) ^ int(R_bit) for R_bit, subkey_bit in zip(p_array_2, subkey_second_half)]

    complete = [first_half, second_half]

    first_half_binary = merge_indices(complete[0])
    second_half_binary = merge_indices(complete[1])

    first_half_decimal_row_column = tuple(binary_to_decimal(pair) for pair in first_half_binary)
    second_half_decimal_row_column = tuple(binary_to_decimal(pair) for pair in second_half_binary)

    bits_from_s0 = s_box_lookup(S0, first_half_decimal_row_column[0], first_half_decimal_row_column[1])
    bits_from_s1 = s_box_lookup(S1, second_half_decimal_row_column[0], second_half_decimal_row_column[1])

    result = [bits_from_s0, bits_from_s1]
    result2 = result[0] + result[1]

    p4_result = p4_function(result2)
    p4_result_2 = p4_result[0] + p4_result[1] + p4_result[2] + p4_result[3]

    return p4_result_2


# Function to perform F-function with subkey for encryption/decryption
def fK(L, R, subkey):
    F_output = F_R_subkey(R, subkey)

    L_xor_F = [int(x) ^ int(y) for x, y in zip(L, F_output)]

    return L_xor_F, R


# Function for S-DES encryption
def encryption(plaintext):
    ip_result = ip_function(plaintext)

    for i in range(len(ip_result)):
        ip_result[i] = int(ip_result[i])

    l_ip_result = ip_result[:4]
    r_ip_result = ip_result[4:]

    result_fk_l, result_fk_r = fK(l_ip_result, r_ip_result, k1)

    switch_result = switch(result_fk_l, result_fk_r)

    fk2_result = fK(switch_result[0], switch_result[1], k2)

    new_list = []
    for row in fk2_result:
        new_list.extend(row)

    ciphertext = ip_inv_function(new_list)

    return ciphertext


# Function for S-DES encryption
def decryption(ciphertext):
    ip_result = ip_function(ciphertext)

    for i in range(len(ip_result)):
        ip_result[i] = int(ip_result[i])

    l_ip_result = ip_result[:4]
    r_ip_result = ip_result[4:]

    result_fk_l, result_fk_r = fK(l_ip_result, r_ip_result, k2)

    switch_result = switch(result_fk_l, result_fk_r)

    fk2_result = fK(switch_result[0], switch_result[1], k1)

    new_list = []
    for row in fk2_result:
        new_list.extend(row)

    plaintext = ip_inv_function(new_list)

    return plaintext


# Perform S-DES encryption and decryption testing
for i in range(1, 11):
    print("\nPlaintext", i, ":", eval("plaintext" + str(i)))
    print("Ciphertext", i, ":", encryption(eval("plaintext" + str(i))))
    print("Plaintext", i, ":", decryption(encryption(eval("plaintext" + str(i)))))

