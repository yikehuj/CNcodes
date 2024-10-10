from traceback import print_stack

from Tools.scripts.summarize_stats import print_title
from pyshamir import split, combine
from setuptools.wheel import unpack

from benchmarks.benchmark_bit_packing import bit_pack_1
from src.dilithium_py.dilithium import *
from src.dilithium_py.modules import *
from src.dilithium_py.modules.modules import MatrixDilithium, ModuleDilithium


class unpack:
    def __init__(self):
        self.M = ModuleDilithium()
        self.l = 4
        self.k = 4
        self.eta = 2

    def unpack_s1(self, bytes_s):
        bytes_s = self.M.bit_unpack_s(bytes_s, self.l, 1, self.eta)
        return bytes_s

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
# p_s2 = split(s2, num_shares, threshold)
# p_s1 = split(s1, num_shares, threshold)
# p_sk = Dilithium2._pack_sk(rho, K, tr, p_s1, p_s2, t0)
# p_sign = Dilithium2.sign(p_sk[0], msg)
# print(p_sign)
# p_s1 = split(s1, num_shares, threshold)

# p_s = MatrixDilithium.bit_pack_s(s1, 2)
# print('bit_s :\n', p_s, '\n')
un_pack = unpack()
# p_s_unpack = un_pack.unpack_s1(p_s)
# print('matrix_s:\n', p_s_unpack,'\n')
# print(s1)
pack_s1 = MatrixDilithium.bit_pack_s(s1, 2)
part_pack_s1 = split(pack_s1, num_shares, threshold)
unpack_s1 = un_pack.unpack_s1(part_pack_s1[0])
print(unpack_s1)
