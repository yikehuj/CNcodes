import os
from src.dilithium_py.modules.modules import ModuleDilithium

try:
    from xoflib import shake256
except ImportError:
    from src.dilithium_py.shake.shake_wrapper import shake256


class Dilithium:
    def __init__(self, parameter_set):
        self.d = parameter_set["d"]
        self.k = parameter_set["k"]
        self.l = parameter_set["l"]
        self.eta = parameter_set["eta"]
        self.tau = parameter_set["tau"]
        self.omega = parameter_set["omega"]
        self.gamma_1 = parameter_set["gamma_1"]
        self.gamma_2 = parameter_set["gamma_2"]
        self.beta = self.tau * self.eta

        self.M = ModuleDilithium()
        self.R = self.M.ring

        # 默认使用系统随机性，实现确定性随机性
        # 使用方法“set_drbg_seed()”
        self.random_bytes = os.urandom

    def set_drbg_seed(self, seed):
        """
        将熵源更改为 DRBG，并用提供的值作为其种子。

        设置种子会将熵源从：func:`os.urandom()`切换到 AES256 CTR DRBG。

        用于 Kyber 的确定性版本以及测试与 KAT 向量的对齐

        Note:
          目前需要 pycryptodome 来实现 AES。
        """
        try:
            from src.dilithium_py.drbg.aes256_ctr_drbg import AES256_CTR_DRBG

            self._drbg = AES256_CTR_DRBG(seed)
            self.random_bytes = self._drbg.random_bytes
        except ImportError as e:  # pragma: no cover
            print(f"Error importing AES from pycryptodome: {e = }")
            raise Warning(
                "由于缺少依赖项，无法设置 DRBG 种子，请尝试安装要求：pip -r install requirements"
            )

    """
    H() 在代码中的几个地方使用 Shake256 将数据哈希为 32 和 64 字节
    """

    @staticmethod
    def _h(input_bytes, length):
        """
        是利用 SHAKE-256 算法生成一个基于输入 input_bytes 的随机或伪随机字节序列，并返回该序列的前 length 个字节
        """
        return shake256(input_bytes).read(length)

    def _expand_matrix_from_seed(self, rho):
        """
        辅助函数从种子“rho”生成大小为k x l 的元素。
        """
        A_data = [[0 for _ in range(self.l)] for _ in range(self.k)]
        for i in range(self.k):
            for j in range(self.l):
                A_data[i][j] = self.R.rejection_sample_ntt_poly(rho, i, j)
        return self.M(A_data)

    def _expand_vector_from_seed(self, rho_prime):
        s1_elements = [
            self.R.rejection_bounded_poly(rho_prime, i, self.eta) for i in range(self.l)
        ]
        s2_elements = [
            self.R.rejection_bounded_poly(rho_prime, i, self.eta)
            for i in range(self.l, self.l + self.k)
        ]

        s1 = self.M.vector(s1_elements)
        s2 = self.M.vector(s2_elements)
        return s1, s2

    def _expand_mask_vector(self, rho_prime, kappa):
        elements = [
            self.R.sample_mask_polynomial(rho_prime, i, kappa, self.gamma_1)
            for i in range(self.l)
        ]
        return self.M.vector(elements)

    @staticmethod
    def _pack_pk(rho, t1):
        return rho + t1.bit_pack_t1()

    def _pack_sk(self, rho, K, tr, s1, s2, t0):
        s1_bytes = s1.bit_pack_s(self.eta)
        s2_bytes = s2.bit_pack_s(self.eta)
        t0_bytes = t0.bit_pack_t0()
        return rho + K + tr + s1_bytes + s2_bytes + t0_bytes

    def _pack_h(self, h):
        non_zero_positions = [
            [i for i, c in enumerate(poly.coeffs) if c == 1]
            for row in h._data
            for poly in row
        ]
        packed = []
        offsets = []
        for positions in non_zero_positions:
            packed.extend(positions)
            offsets.append(len(packed))

        padding_len = self.omega - offsets[-1]
        packed.extend([0 for _ in range(padding_len)])
        return bytes(packed + offsets)

    def _pack_sig(self, c_tilde, z, h):
        return c_tilde + z.bit_pack_z(self.gamma_1) + self._pack_h(h)

    def _unpack_pk(self, pk_bytes):
        rho, t1_bytes = pk_bytes[:32], pk_bytes[32:]
        t1 = self.M.bit_unpack_t1(t1_bytes, self.k, 1)
        return rho, t1

    def _unpack_sk(self, sk_bytes):
        if self.eta == 2:
            s_bytes = 96
        else:
            s_bytes = 128
        s1_len = s_bytes * self.l
        s2_len = s_bytes * self.k
        t0_len = 416 * self.k
        # if len(sk_bytes) != 3 * 32 + s1_len + s2_len + t0_len:
        #     raise ValueError("SK packed bytes is of the wrong length")

        # Split bytes between seeds and vectors
        sk_seed_bytes, sk_vec_bytes = sk_bytes[:96], sk_bytes[96:]

        # Unpack seed bytes
        rho, K, tr = (
            sk_seed_bytes[:32],
            sk_seed_bytes[32:64],
            sk_seed_bytes[64:96],
        )

        # Unpack vector bytes
        s1_bytes = sk_vec_bytes[:s1_len]
        s2_bytes = sk_vec_bytes[s1_len : s1_len + s2_len]
        t0_bytes = sk_vec_bytes[-t0_len:]

        # Unpack bytes to vectors
        s1 = self.M.bit_unpack_s(s1_bytes, self.l, 1, self.eta)
        s2 = self.M.bit_unpack_s(s2_bytes, self.k, 1, self.eta)
        t0 = self.M.bit_unpack_t0(t0_bytes, self.k, 1)

        return rho, K, tr, s1, s2, t0

    def _unpack_h(self, h_bytes):
        offsets = [0] + list(h_bytes[-self.k :])
        non_zero_positions = [
            list(h_bytes[offsets[i] : offsets[i + 1]]) for i in range(self.k)
        ]

        matrix = []
        for poly_non_zero in non_zero_positions:
            coeffs = [0 for _ in range(256)]
            for non_zero in poly_non_zero:
                coeffs[non_zero] = 1
            matrix.append([self.R(coeffs)])
        return self.M(matrix)

    def _unpack_sig(self, sig_bytes):
        c_tilde = sig_bytes[:32]
        z_bytes = sig_bytes[32 : -(self.k + self.omega)]
        h_bytes = sig_bytes[-(self.k + self.omega) :]

        z = self.M.bit_unpack_z(z_bytes, self.l, 1, self.gamma_1)
        h = self._unpack_h(h_bytes)
        return c_tilde, z, h

    def keygen(self):
        """
        Generates a public-private keyair
        """
        # Random seed
        zeta = self.random_bytes(32)

        # Expand with an XOF (SHAKE256)
        seed_bytes = self._h(zeta, 128)

        # Split bytes into suitable chunks
        rho, rho_prime, K = seed_bytes[:32], seed_bytes[32:96], seed_bytes[96:]

        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Generate the error vectors s1 ∈ R^l, s2 ∈ R^k
        s1, s2 = self._expand_vector_from_seed(rho_prime)
        s1_hat = s1.to_ntt()

        # Matrix multiplication
        t = (A_hat @ s1_hat).from_ntt() + s2

        t1, t0 = t.power_2_round(self.d)

        # Pack up the bytes
        pk = self._pack_pk(rho, t1)
        tr = self._h(pk, 32)

        sk = self._pack_sk(rho, K, tr, s1, s2, t0)
        return pk, sk
# 公钥pk中，ρ是哈希函数H的输入，t1是签名的一部分
# 私钥sk中，ρ是哈希函数H的输入，K是私钥，tr是另一个哈希函数H’的输出，s1和s2是随机向量，t0是签名的一部分。
# 在后续的签名和验证过程中，私钥将用于生成签名，公钥将用于验证签名。



    def sign(self, sk_bytes, m):
        """
        Generates a signature for a message m from a byte-encoded private key
        """
        # unpack the secret key
        rho, K, tr, s1, s2, t0 = self._unpack_sk(sk_bytes)

        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Set seeds and nonce (kappa)
        mu = self._h(tr + m, 64)
        kappa = 0
        rho_prime = self._h(K + mu, 64)

        # Precompute NTT representation
        s1 = s1.to_ntt()
        s2 = s2.to_ntt()
        t0 = t0.to_ntt()

        alpha = self.gamma_2 << 1
        while True:
            y = self._expand_mask_vector(rho_prime, kappa)
            y_hat = y.to_ntt()

            # increment the nonce
            kappa += self.l

            w = (A_hat @ y_hat).from_ntt()

            # Extract out both the high and low bits
            w1, w0 = w.decompose(alpha)

            # Create challenge polynomial
            w1_bytes = w1.bit_pack_w(self.gamma_2)
            c_tilde = self._h(mu + w1_bytes, 32)
            c = self.R.sample_in_ball(c_tilde, self.tau)

            # Store c in NTT form
            c = c.to_ntt()

            z = y + (s1.scale(c)).from_ntt()
            if z.check_norm_bound(self.gamma_1 - self.beta):
                continue

            w0_minus_cs2 = w0 - s2.scale(c).from_ntt()
            if w0_minus_cs2.check_norm_bound(self.gamma_2 - self.beta):
                continue

            c_t0 = t0.scale(c).from_ntt()
            if c_t0.check_norm_bound(self.gamma_2):
                continue

            w0_minus_cs2_plus_ct0 = w0_minus_cs2 + c_t0

            h = w0_minus_cs2_plus_ct0.make_hint_optimised(w1, alpha)
            if h.sum_hint() > self.omega:
                continue

            return self._pack_sig(c_tilde, z, h)
        # 这行代码是Dilithium签名算法的输出，返回一个签名 σ其中，c是签名中的哈希值的NTT表示，
        # z是签名中的向量y和一些其他值的线性组合，h是生成签名的过程中计算出的一些用于验证签名的值的哈希值的NTT表示。


    def verify(self, pk_bytes, m, sig_bytes):
        """
        Verifies a signature for a message m from a byte encoded public key and
        signature
        """
        rho, t1 = self._unpack_pk(pk_bytes)
        c_tilde, z, h = self._unpack_sig(sig_bytes)

        if h.sum_hint() > self.omega:
            return False

        if z.check_norm_bound(self.gamma_1 - self.beta):
            return False

        A_hat = self._expand_matrix_from_seed(rho)

        tr = self._h(pk_bytes, 32)
        mu = self._h(tr + m, 64)
        c = self.R.sample_in_ball(c_tilde, self.tau)

        # Convert to NTT for computation
        c = c.to_ntt()
        z = z.to_ntt()

        t1 = t1.scale(1 << self.d)
        t1 = t1.to_ntt()

        Az_minus_ct1 = (A_hat @ z) - t1.scale(c)
        Az_minus_ct1 = Az_minus_ct1.from_ntt()

        w_prime = h.use_hint(Az_minus_ct1, 2 * self.gamma_2)
        w_prime_bytes = w_prime.bit_pack_w(self.gamma_2)

        return c_tilde == self._h(mu + w_prime_bytes, 32)
