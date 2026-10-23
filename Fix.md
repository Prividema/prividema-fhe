# Code review findings

Start with `src/core/polynomial/bivariate_polynomial.c:347`. It's the bug most likely to hit you quietly, and it was confirmed by running the code. The project was built and a small probe program was run against it. The first three bugs below are confirmed by that run; the rest come from reading the code.

## 1. Confirmed bugs (checked by running the code)

1. **`biv_decomp_internal` leaves one limb uninitialized** (`bivariate_polynomial.c:347-351`). The final zeroing loop starts at `last_l - nb_limbs_to_write - 2`, so limb `last_l - nb_limbs_to_write - 1` is never written. Probe: `univ_rnx_to_biv` on 2^-40 into a buffer pre-filled with 12345 gave limbs `12345 0 131072 0 0`. Any reused buffer (encryption temp buffers, for example) carries garbage into the result.
2. **HEInt is wrong for large `t`** (`heint.c:49-61`, `103-115`). The sums are added in `uint32_t`, which wraps mod 2^32, not mod t. Probe: N=1024 with t=1073750017 (which is < 2^31) gave 1024/1024 wrong slots. `heint.h:24` says any t ≤ 2^31 works. The tests only pass because 1024·4085761 is just under 2^32.
3. **`glwe_params_stdev` scales NOISE_NORMAL twice** (`glwe_params.c:40` vs `:76`). The constructor stores the raw sigma, and stdev then multiplies it by 2^-Kℓa again. Probe: sigma 2^-60 came back as 2^-155.

## 2. Memory and logic bugs (from reading)

1. **Double free on error paths**: `ggsw_arithmetic.c:256` + `264-267` and `glwegadget_arithmetic.c:533` + `541-544`. Each element is freed inside the loop but not set to NULL, then the cleanup block frees all of them again.
2. **`printf_poly_biv` reads out of bounds** (`bivariate_polynomial.c:784`). It reads `ptr[(i - 1) * stride]` with `i` starting at 0, which reads before the buffer.
3. **GGSW code uses `k` where the docs say `k_tilde`** (`ggsw_ciphertext.c:97-102`, `ggsw_arithmetic.c:178`). The header (`ggsw_ciphertext.h:63`) and `ggsw_glwe_view` use `k_tilde`. With `k_tilde < k`, which the docs allow, rows land in the wrong place.
4. **`hebits_cmux_tree` with `inp_cols == 1`** (`hebits.c:54,72`). The loop runs 0 times, so `res` is never written and the function still returns 0. With `delete_src`, the input also leaks.
5. **`glwe_trace_expand` with `res_size == 1`** (`glwegadget_arithmetic.c:333`). It writes `results[1]`, which is out of bounds. It also scales by 2, while the header formula `m = 2^ceil(log2(res_size))` gives 1.

## 3. Comments that contradict the code

1. **`univ_tnx_to_biv` rounding**: the header (`bivariate_polynomial.h:189`) says the dropped LSB is "rounded down, not towards 0". The code (`.c:495-497`) takes the absolute value and then shifts, which rounds towards 0.
2. **`hefp_encode`/`hefp_decode` scale sign**: the docs say Δ = 2^-scale_bits (`hefp.h:22`). The code passes `-scale_bits` as a right-shift count, so Δ = 2^+scale_bits. The tests use scale_bits ≤ 0, which matches the code, not the docs.
3. **`glwegadget_key_collection_put_key`** returns `(void*)-1` when `pos` is out of range (`glwegadget_key.c:104`). The header only mentions "previous occupant or NULL".
4. **NOISE_NORMAL is accepted but can't be used**. `new_glwe_params` accepts it, but `biv_add_noise` rejects it (`bivariate_polynomial.c:170`) and `rand_normal` always fails (`rng.c:23`). `bench_glwe_encrypt.cpp:82` uses it, so its encryptions fail.
5. **`glwegadget_automorphism`**: the comment says the precision comes from "the input b" (`:172-175`). The code uses the *result's* `l_b`.

## 4. Inconsistencies that can hide failures

1. **Errors reported as success**:
   - These always return 0, even when the work failed: `glwe_dft_mult_const` (`glwe_arithmetic.c:101`), `glwe_mult_const` (`:78`), `spqlios_vec_znx_idft` (`spqlios_alias.c:46-54`, even when its malloc fails), and `hefp_*_internal`.
   - These ignore the return value of the call they make: `heint_encode` and `heint_decode` (`heint.c:161,207`) and the trace-expand functions.
2. **`ggsw_secret_encrypt` can report success on a parameter error** (`ggsw_ciphertext.c:73`). It returns `log_perror(...)`, which is 0 when LOG_ERROR is masked. The message would also print an unrelated `strerror(errno)`.
3. **Some objects skip the backend allocators**: `new_biv_dft` (`calloc`), `new_ggsw_prep` (`calloc`, no `module`), and `ggsw_external_product_unprepared` (it frees a buffer with `delete_biv_dft`). `GLWEGadgetCiphertextPrep` and `new_univ_dft` go through the backend. This is the same class of problem as commit a984fd8 (mismatched frees).

## 5. Cosmetic
  - Some messages still use old names: `glwe_secret_masking_ggsw_lib`, `new_biv_poly`.

## 6. Unit tests to implement

The current tests miss these bugs because they only use clean buffers, small parameters and the default code paths. Each test below should fail on the current code and pass once the fix is in. Write each one before its fix so you can see it go from red to green. All of them fit into the existing Criterion files under `tests/unit/`.

1. **Stale limb in `univ_rnx_to_biv`** (`core/glwe/test_bivariate_polynomial.c`): fill the output `PolyBiv` with a non-zero pattern, convert values like 2^-40 and 0.0, then check that every limb matches a conversion into a zeroed buffer.
2. **HEInt with a large `t`** (`schemes/test_heint_encdec.c`): run encode then decode with N=1024 and t=1073750017 (prime, t ≡ 1 mod 2N), plus one t close to 2^31. Check that all slots come back unchanged.
3. **`glwe_params_stdev` for NOISE_NORMAL** (`core/glwe/test_glwe_ct_params.c`): create params with sigma = 2^-60 and check that the returned stdev is 2^-60.
4. **Error-path frees** (`core/ggsw/test_trace_expand.c`): make the prepared trace expansions fail partway, for example with a missing automorphism key. Check that they return -1, and run under AddressSanitizer to catch the double free.
5. **Edge sizes** (`schemes/test_hebits.c`, `core/ggsw/test_trace_expand.c`): `hebits_cmux_tree` with `inp_cols == 1` should copy `src[0]` into `res`. `glwe_trace_expand` with `res_size == 1` should not write `results[1]`.
6. **`k_tilde < k`** (`core/ggsw/test_ggsw_encrypt.c`): encrypt a GGSW with k=2 and k_tilde=1. Decrypt each row through `ggsw_glwe_view` and check that it matches the documented layout.
7. **Documented behaviour** (`core/glwe/test_bivariate_polynomial.c`, `schemes/test_hefp_encoding.c`, `core/ggsw/test_automorphisms.c`):
   - `univ_tnx_to_biv` on -1 should round the way the header says.
   - `hefp_encode` with scale_bits = -9: decide whether Δ should be 2^-9 (what the code does) or 2^9 (what the header says), then assert that value.
   - `glwegadget_key_collection_put_key` with an out-of-range `pos` should return the documented value.

   Fix either the code or the comment, then keep the test.

Fix time: about 15–30 minutes each for the groups 1–2 bugs. The comment fixes take under an hour in total.

Next: open `bivariate_polynomial.c:347` and change `nb_limbs_to_write + 1` to `nb_limbs_to_write` in the `min_i_2` line. Then rerun the probe to check that limb 0 prints 0.

# TODO
Rename/Refactor ZnxVec.