from pyshamir import split, combine
from src.dilithium_py.dilithium import Dilithium2
pk, sk = Dilithium2.keygen()
# 示例用法
rho, K, tr, s1, s2, t0 = Dilithium2._unpack_sk(sk)
print('secret:\n', sk)
msg = b'im winner'
threshold = 4  # 门限值
num_shares = 5  # 份额数量
# secret = secret.to_bytes(32)
sign = Dilithium2.sign(sk, msg)
parts = split(sign, num_shares, threshold)
for i in range(len(parts)):
    print('parts[', i + 1, ']:\n', parts[i], '\n')
recomb_secret = combine(parts)
print('recomb_secret:\n', recomb_secret)
print('sig验证:\n', Dilithium2.verify(pk, msg, recomb_secret))

