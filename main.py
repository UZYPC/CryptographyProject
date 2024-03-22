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


