#pragma once
#ifndef YDAR95_GALOIS_FIELD_H
#define YDAR95_GALOIS_FIELD_H
#include <cstdint>
#include <atomic>
/**
 * \brief Math of Galois Field ,base=2 GF(2^BitWide) f(x)=Prim
 * \tparam BitWide	位数 ，	取值[1,12]
 * \tparam Prim		多项式，	取值为质数，且Prim_binary[BitWide]=1; Prim_binary[ >BitWide ]=0 , defaut=0，使用预设值
 */
template<uint32_t BitWide,uint32_t Prim=0>
class GaloisField
{
public:
	GaloisField(const uint32_t Alpha=0)
	:alpha(Alpha)
	{
		static_assert(BitWide > 0 && BitWide <= 12,"BitWide between 1 and 12");
		if (sm_table_is_init == false)
		{
			sm_table_is_init = true;
			gf_init();
		}
	}
	~GaloisField()
	{
		
	}

	GaloisField<BitWide, Prim> neg() const
	{
		return *this;
	}
	GaloisField<BitWide, Prim> inv() const
	{
		return GaloisField<BitWide, Prim>(sm_table_div[1][this->alpha]);
	}

	friend GaloisField<BitWide, Prim> operator+ (
		const GaloisField<BitWide, Prim>& lhs, const GaloisField<BitWide, Prim>& rhs)
	{
		return  GaloisField<BitWide, Prim>(gf_add(lhs.alpha, rhs.alpha));
	}
	friend GaloisField<BitWide, Prim> operator- (
		const GaloisField<BitWide, Prim>& lhs, const GaloisField<BitWide, Prim>& rhs)
	{
		return  GaloisField<BitWide, Prim>(gf_sub(lhs.alpha, rhs.alpha));
	}
	friend GaloisField<BitWide, Prim> operator* (
		const GaloisField<BitWide, Prim>& lhs, const GaloisField<BitWide, Prim>& rhs)
	{
		return  GaloisField<BitWide, Prim>(gf_mul(lhs.alpha, rhs.alpha));
	}
	friend GaloisField<BitWide, Prim> operator/ (
		const GaloisField<BitWide, Prim>& lhs, const GaloisField<BitWide, Prim>& rhs)
	{
		return  GaloisField<BitWide, Prim>(gf_div(lhs.alpha, rhs.alpha));
	}
public:
	static const uint32_t MAX_VAL;// = (1 << BitWide) - 1;
	static const uint32_t BIT_WIDE;
public:
	uint32_t alpha;

protected:
	static uint32_t gf_mul(uint32_t a, uint32_t b)
	{
		return sm_table_mul[a][b];
	}
	static uint32_t gf_div(uint32_t a, uint32_t b)
	{
		return sm_table_div[a][b];
	}
	static uint32_t gf_add(uint32_t a, uint32_t b)
	{
		return a^b;
	}
	static uint32_t gf_sub(uint32_t a, uint32_t b)
	{
		return a^b;
	}

	static uint32_t gf_no_table_mul(uint32_t a, uint32_t b)
	{
		if (0 == a || 0 == b)
			return 0;
		return sm_table_alpha[(sm_table_index[a] + sm_table_index[b]) % (MAX_VAL)];
	}
	static uint32_t gf_no_table_div(uint32_t a, uint32_t b)
	{
		if (0 == a || 0 == b)
			return 0;
		return sm_table_alpha[(sm_table_index[a] - sm_table_index[b] + MAX_VAL) % (MAX_VAL)];
	}
	static void     gf_init()
	{
		const uint32_t prim_poly[13] =
		{
			/*	0 */	0x00000000,
			/*  1 */    0x00000001,
			/*  2 */    0x00000007,
			/*  3 */    0x0000000b,
			/*  4 */    0x00000013,
			/*  5 */    0x00000025,
			/*  6 */    0x00000043,
			/*  7 */    0x00000089,
			/*  8 */    0x00000187,
			/*  9 */    0x00000211,
			/* 10 */    0x00000409,
			/* 11 */    0x00000805,
			/* 12 */    0x00001053,
		};
		uint32_t GFprim = prim_poly[BitWide];
		if (Prim != 0)
			GFprim = Prim;

		const int      GFsize = 1 << BitWide;
		sm_table_alpha[0] = 1;
		sm_table_index[0] = -1;// => MAX
		for (int i = 1; i<GFsize; i++)
		{
			sm_table_alpha[i] = sm_table_alpha[i - 1] << 1;
			if (sm_table_alpha[i] >= GFsize)
			{
				sm_table_alpha[i] ^= GFprim;
			}
			sm_table_index[sm_table_alpha[i]] = i;
		}
		sm_table_index[1] = 0;

		for (int i = 0; i < GFsize; i++) {
			for (int j = 0; j < GFsize; j++)
			{
				sm_table_mul[i][j] = gf_no_table_mul(i, j);
				sm_table_div[i][j] = gf_no_table_div(i, j);
			}
		}
	}

protected:
	static uint32_t  sm_table_alpha[1<<BitWide];
	static uint32_t  sm_table_index[1 << BitWide];
	static uint32_t  sm_table_mul[1 << BitWide][1 << BitWide];
	static uint32_t  sm_table_div[1 << BitWide][1 << BitWide];
	static std::atomic<bool> sm_table_is_init;
};



template<uint32_t BitWide, uint32_t Prim>
const uint32_t  GaloisField<BitWide, Prim>::BIT_WIDE = BitWide;

template<uint32_t BitWide, uint32_t Prim>
const uint32_t  GaloisField<BitWide, Prim>::MAX_VAL = (1 << BitWide) - 1;

template<uint32_t BitWide, uint32_t Prim>
uint32_t  GaloisField<BitWide, Prim>::sm_table_alpha[1 << BitWide];

template<uint32_t BitWide, uint32_t Prim>
uint32_t  GaloisField<BitWide, Prim>::sm_table_index[1 << BitWide];

template<uint32_t BitWide, uint32_t Prim>
uint32_t  GaloisField<BitWide, Prim>::sm_table_mul[1 << BitWide][1 << BitWide];

template<uint32_t BitWide, uint32_t Prim>
uint32_t  GaloisField<BitWide, Prim>::sm_table_div[1 << BitWide][1 << BitWide];

template<uint32_t BitWide, uint32_t Prim>
std::atomic<bool> GaloisField<BitWide, Prim>::sm_table_is_init;

#endif //YDAR95_GALOIS_FIELD_H