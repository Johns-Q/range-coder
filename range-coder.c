///
///	@file range-coder.c	@brief RangeCoder: Range encoding
///
///	Copyright (c) 2015 by Johns.	 All Rights Reserved.
///
///	Contributor(s):
///
///	License: AGPLv3
///
///	This program is free software: you can redistribute it and/or modify
///	it under the terms of the GNU Affero General Public License as
///	published by the Free Software Foundation, either version 3 of the
///	License.
///
///	This program is distributed in the hope that it will be useful,
///	but WITHOUT ANY WARRANTY; without even the implied warranty of
///	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
///	GNU Affero General Public License for more details.
///
///	$Id$
//////////////////////////////////////////////////////////////////////////////

///
///	@defgroup RangeCoder	Range encoding
///
///	https://en.wikipedia.org/wiki/Range_encoding
///
/// @{

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

/**
**	RangeCoder typedef.
*/
typedef struct _range_coder_ RangeCoder;

unsigned LoopsMax;

#if !defined(RANGE_CODER_32) && !defined(RANGE_CODER_64)

#define RANGE_CODER_64			///< full 64-bit version

#endif

#if !defined(RANGE_CODER_32) && !defined(RANGE_CODER_64)

#define RANGE_CODER_32			///< full 32-bit version

#endif

///
///	@addtogroup RangeCoder
///
///	- 64-bit
///
///	7	 6	  5	   4	    3	     2	      1	       0
///	bbbbbbbb|pppppppp|puuuuuuu|uuuuuuuu|uuuuuuuu|uuuuuuuu|uuuuuuuu|uuuuuuuu
///			   ^- max usable sum
///	       ^ RangeCoderRangeBot => RangeCoderShift
///	^ RangeCoderRangeTop => RangeCoderCodeBits
enum
{
#ifdef RANGE_CODER_64
    RangeCoderCodeBits = 64,
    // too big: RangeCoderTotSum = UINT64_C(1) << (64 - 8 - 9),
    RangeCoderTotSum = UINT64_C(1) << 31,
    RangeCoderShift = RangeCoderCodeBits - 8,
    RangeCoderRangeBot = UINT64_C(1) << RangeCoderShift,
    RangeCoderRangeTop = UINT64_C(-1),
#endif
#ifdef RANGE_CODER_32
    RangeCoderCodeBits = 32,
    RangeCoderTotSum = UINT32_C(1) << (32 - 8 - 9),
    RangeCoderShift = RangeCoderCodeBits - 8,
    RangeCoderRangeBot = UINT32_C(1) << RangeCoderShift,
    RangeCoderRangeTop = UINT32_C(-1),
#endif
};

#ifdef RANGE_CODER_64
typedef uint64_t rc_uint_t;		///< range-coder range integer
#endif
#ifdef RANGE_CODER_32
typedef uint32_t rc_uint_t;		///< range-coder range integer
#endif

/**
**	RangeCoder encoder and decoder context.
*/
struct _range_coder_
{
    rc_uint_t Low;			///< start of range
    rc_uint_t Range;			///< total range
    rc_uint_t Save;			///< decoder save result range / totsum

    unsigned Buffer;			///< 1 byte encoder buffer
    unsigned Overflow;			///< encoder overflow counter

    uint8_t *Data;			///< in/out data buffer
    unsigned Size;			///< data buffer size
    unsigned Used;			///< used area of data buffer
    unsigned Index;			///< index into data buffer
    ///
    /// Function to stream data buffer in/output
    ///
    unsigned (*Stream) (void *, uint8_t **, unsigned *);
    void *User;				///< user data for stream function
};

#if GCC_VERSION >= 50100
///
///	Add with overflow check, added in future GCC >=5.1 versions.
///
///	@TODO Can be used for carry detection
///
#define check_add_overflow(a, b, d) ({		\
	typeof(a) __a = (a);			\
	typeof(b) __b = (b);			\
	typeof(d) __d = (d);			\
	(void) (&__a == &__b);			\
	(void) (&__a == __d);			\
	__builtin_add_overflow(__a, __b, __d);	\
})
#else
#endif

/**
**	Stream reader/writer which asserts, if called.
**
**	@param user	opaque user data pointer
**	@param data	address of buffer pointer
**	@param size	address of buffer size
**
**	@returns index to continue in buffer
*/
unsigned RangeCoderStreamAssert( __attribute__ ((unused))
    void *user, __attribute__ ((unused)) uint8_t ** data, unsigned *size)
{
    assert(!*size);
    return *size - 1;
}

/**
**	No operation stream writer.
**
**	@param user	opaque user data pointer
**	@param data	address of buffer pointer
**	@param size	address of buffer size
**
**	@returns index to end of buffer.
*/
unsigned RangeCoderWriteNoop( __attribute__ ((unused))
    void *user, __attribute__ ((unused)) uint8_t ** data, unsigned *size)
{
    return *size - 1;
}

/**
**	Realloc buffer stream writer.
**
**	@param user	opaque user data pointer
**	@param data	address of buffer pointer
**	@param size	address of buffer size
**
**	realloc buffer, return *size - 1, if failed.
*/
unsigned RangeCoderWriteRealloc( __attribute__ ((unused))
    void *user, uint8_t ** data, unsigned *size)
{
    void *tmp;
    unsigned n;

    n = *size;
    if ((tmp = realloc(*data, n + 8192))) {
	*data = tmp;
	*size = n + 4096;
	return n;
    }
    return n - 1;
}

/**
**	Stdio buffer stream writer.
**
**	@param user	opaque user data pointer
**	@param data	address of buffer pointer
**	@param size	address of buffer size
**
**	@returns *size-1, if failed.
*/
unsigned RangeCoderWriteStdio( __attribute__ ((unused))
    void *user, uint8_t ** data, unsigned *size)
{
    size_t n;

    n = fwrite(*data, 1, *size, (FILE *) user);
    if (n == *size) {
	return *size - n;
    }
    return *size - 1;
}

/**
**	Stdio buffer stream reader.
**
**	@param user	opaque user data pointer
**	@param data	address of buffer pointer
**	@param size	address of buffer size
**
**	@returns number of bytes read into buffer.
*/
unsigned RangeCoderReadStdio( __attribute__ ((unused))
    void *user, uint8_t ** data, unsigned *size)
{
    size_t n;

    n = fread(*data, 1, *size, (FILE *) user);

    return n;
}

/**
**	Initialize range-coder for encoding.
**
**	@param rc		RangeCoder context (NULL allocate new)
**	@param first_byte	byte, which is always written unencoded
**	@param data		encode buffer
**	@param size		size of encode buffer
**	@param stream		function to write buffer, called if buffer full
**	@param user		argument for stream write function
**
**	@retval NULL	out of memory
**	@returns pointer to RangeCoder context.
*/
RangeCoder *RangeCoderEncodeInit(RangeCoder * rc, unsigned first_byte,
    uint8_t * data, unsigned size, unsigned (*stream) (void *, uint8_t **,
	unsigned *), void *user)
{
    if (!rc) {
	if (!(rc = malloc(sizeof(*rc)))) {
	    return NULL;		// out of memory
	}
    }
    // prepare output stream
    rc->Data = data;
    rc->Size = size;
    rc->Stream = stream;
    rc->User = user;

    if (!rc->Data && !rc->Stream) {	// no buffer and no writer
	rc->Data = (uint8_t *) & rc->Buffer;
	rc->Stream = RangeCoderWriteNoop;
	rc->Size = 1;
    }
    if (!rc->Stream) {			// use dummy to avoid NULL check
	rc->Stream = RangeCoderStreamAssert;
    }
    rc->Index = 0;

    rc->Low = 0U;			// setup encoder
    rc->Buffer = first_byte;
    rc->Overflow = 1U;
    rc->Range = RangeCoderRangeTop;

    return rc;
}

/**
**	Write byte to output stream.
**
**	@param rc	RangeCoder context
**	@param byte	encoded byte to write to output stream
*/
static inline void EncodeWriteByte(RangeCoder * rc, unsigned byte)
{
    assert(rc->Data);

    if (rc->Index >= rc->Size) {
	rc->Index = rc->Stream(rc->User, &rc->Data, &rc->Size);
	// full sets index = size - 1
    }
    assert(rc->Index < rc->Size);
    rc->Data[rc->Index++] = byte;
}

/**
**	Renormalize the encoder.
**
**	@param rc	RangeCoder context
*/
static inline void EncodeRenormalize(RangeCoder * rc)
{
    unsigned loops;

    assert(rc->Range <= RangeCoderRangeTop);

    // range has become too small
    loops = 0;
    // FIXME: we can unroll loop, max depending on RangeCoderTotSum
    while (rc->Range < RangeCoderRangeBot) {
	assert(rc->Range);

	if (rc->Low < (UINT64_C(0xFF) << RangeCoderShift)) {	// no overflow
	    EncodeWriteByte(rc, rc->Buffer);
	    while (--rc->Overflow) {	// output overflow queue
		EncodeWriteByte(rc, 0xFF);
	    }
	    rc->Buffer = rc->Low >> RangeCoderShift;
	}

	rc->Overflow++;
	rc->Range <<= 8;
	rc->Low <<= 8;
	if (++loops > LoopsMax) {
	    LoopsMax = loops;
	}
    }

    assert(rc->Range <= RangeCoderRangeTop);
}

/**
**	Encode symbol.
**
**	@param rc	RangeCoder context
**	@param lowsum	sum of symbols below encoded symbol
**	@param symsum	sum of encoded symbol
**	@param totsum	sum of all symbols
*/
void RangeCoderEncode(RangeCoder * rc, unsigned lowsum, unsigned symsum,
    unsigned totsum)
{
    rc_uint_t r;
    rc_uint_t l;
    unsigned carry;

    assert(symsum != 0);
    assert(lowsum + symsum <= totsum);
    assert(totsum <= RangeCoderTotSum);

    r = rc->Range / totsum;
    l = r * lowsum;
    carry = (l + rc->Low) < rc->Low;
    rc->Low += l;
    if (lowsum + symsum >= totsum) {	// hit top
	rc->Range -= l;
    } else {
	rc->Range = r * symsum;
    }
    //printf("%08jx\n", (uintmax_t)rc->Range);

    if (carry) {			// output carry early
	rc->Buffer++;
	while (--rc->Overflow) {	// output overflow queue
	    EncodeWriteByte(rc, rc->Buffer);
	    rc->Buffer = 0x00;
	}
	++rc->Overflow;
    }

    EncodeRenormalize(rc);
}

/**
**	Flush the bytes still in encoder buffer.
**
**	@param rc		RangeCoder context
**	@param[out] buffer	returns address of data buffer
**
**	@returns number of bytes in data buffer.
*/
unsigned RangeCoderEncodeFlush(RangeCoder * rc, uint8_t ** buffer)
{
    int i;
    unsigned carry;

    carry = (rc->Range - 1 + rc->Low) < rc->Low;
    rc->Low += rc->Range - 1;

    // output buffer and add carry
    EncodeWriteByte(rc, rc->Buffer + carry);
    while (--rc->Overflow) {		// output overflow queue
	EncodeWriteByte(rc, 0xFF + carry);
    }

    // output remaining code bytes
    EncodeWriteByte(rc, rc->Low >> RangeCoderShift);

    // decoder uses this, than we need no eof check
    for (i = 1; i < RangeCoderCodeBits / 8; ++i) {
	EncodeWriteByte(rc, 0x00);
    }

    if (buffer) {			// report final data buffer
	*buffer = rc->Data;
    }
    return rc->Index;			// report how many in buffer
}

/**
**	Read byte from stream.
**
**	@param rc		RangeCoder context
**
**	@returns next byte from data buffer.
*/
static inline unsigned DecodeReadByte(RangeCoder * rc)
{
    if (rc->Index >= rc->Used) {
	rc->Data[0] = 0x00;		// read 0, if EOF or error
	rc->Used = rc->Stream(rc->User, &rc->Data, &rc->Size);
	rc->Index = 0;
    }
    return rc->Data[rc->Index++];
}

/**
**	Initialize range-coder for decoding.
**
**	@param rc		RangeCoder context (NULL allocate new)
**	@param first_byte	byte, which is always read unencoded
**	@param data		decode buffer
**	@param size		size of decode buffer
**	@param used		bytes used in decode buffer
**	@param stream		function to read buffer, called if buffer empty
**	@param user		argument for stream read function
**
**	@retval NULL	out of memory
**	@returns pointer to RangeCoder context.
*/
RangeCoder *RangeCoderDecodeInit(RangeCoder * rc, unsigned *first_byte,
    uint8_t * data, unsigned size, unsigned used, unsigned (*stream) (void *,
	uint8_t **, unsigned *), void *user)
{
    unsigned byte;
    int i;

    if (!rc) {
	if (!(rc = malloc(sizeof(rc)))) {
	    return NULL;		// out of memory
	}
    }
    // prepare output stream
    rc->Data = data;
    rc->Used = used;
    rc->Size = size;
    rc->Stream = stream;
    rc->User = user;
    if (!rc->Data && !rc->Stream) {	// no buffer or input stream
	rc->Data = (uint8_t *) & rc->Buffer;
	rc->Used = 1;
	rc->Size = 1;
    }
    if (!rc->Stream) {			// dummy stream to avoid null check
	rc->Stream = RangeCoderStreamAssert;
    }
    rc->Index = 0;

    byte = DecodeReadByte(rc);		// unused first byte
    if (first_byte) {
	*first_byte = byte;
    }

    rc->Low = 0;			// fill low
    for (i = 0; i < RangeCoderCodeBits / 8; ++i) {
	rc->Low <<= 8;
	rc->Low += DecodeReadByte(rc);
    }
    rc->Range = RangeCoderRangeTop;

    return rc;
}

/**
**	Renormalize the decoder.
**
**	@param rc		RangeCoder context (NULL allocate new)
*/
static inline void DecodeRenormalize(RangeCoder * rc)
{
    assert(rc->Range <= RangeCoderRangeTop);

    // range has become too small
    while (rc->Range < RangeCoderRangeBot) {
	assert(rc->Range);

	rc->Low <<= 8;
	rc->Low += DecodeReadByte(rc);
	rc->Range <<= 8;
    }

    assert(rc->Range <= RangeCoderRangeTop);
}

/**
**	Get the next symbol sum.
**
**	@param rc	RangeCoder context
**	@param totsum	sum of all symbols
*/
unsigned RangeCoderDecode(RangeCoder * rc, unsigned totsum)
{
    rc_uint_t l;

    assert(totsum <= RangeCoderTotSum);

    rc->Save = rc->Range / totsum;
    l = rc->Low / rc->Save;
    if (l < totsum) {
	return l;
    }
    return totsum - 1;
}

/**
**	Update decoder for symbol.
**
**	@param rc	RangeCoder context
**	@param lowsum	sum of symbols below decoded symbol
**	@param symsum	sum of decoded symbol
**	@param totsum	sum of all symbols
*/
void RangeCoderDecodeUpdate(RangeCoder * rc, unsigned lowsum, unsigned symsum,
    unsigned totsum)
{
    rc_uint_t l;

    assert(symsum != 0);
    assert(lowsum + symsum <= totsum);
    assert(totsum <= RangeCoderTotSum);

    l = rc->Save * lowsum;
    rc->Low -= l;
    if (lowsum + symsum >= totsum) {	// hit top
	rc->Range -= l;
    } else {
	rc->Range = rc->Save * symsum;
    }

    DecodeRenormalize(rc);
}

// ----------------------------------------------------------------------------
///
///	@defgroup bit	 BIT: Binary index tree.
///
///	https://en.wikipedia.org/wiki/Fenwick_tree
///
/// @{

/**
**	Binary indexed tree.
*/
static void BitInit(unsigned *bit_array, int size)
{
    int i;

    for (i = 0; i < size; ++i) {
	bit_array[i] = 0;
    }
}

/**
**	Binary index tree.
*/
static void BitUpdate(unsigned *bit_array, int size, int i, int delta)
{
    assert(i);
    while (i < size) {
	bit_array[i] += delta;
	i += i & -i;
    }
}

/**
**	Divide each original frequency by div.
*/
void BitScaleDown(unsigned *bit_array, int size, int div)
{
    int i;

    assert(div != 0);
    for (i = 1; i < size; i++) {
	bit_array[i] = bit_array[i] / div;
    }
}

/**
**	Get the sum from range (1 ... i) from Binary index tree.
*/
static int BitGetSum(const unsigned *bit_array, int i)
{
    int sum;

    sum = 0;
    while (i > 0) {
	sum += bit_array[i];
	i -= i & -i;
    }
    return sum;
}

/**
**	Binary index tree.
*/
static int BitGet(const unsigned *bit_array, int i)
{
    int sum;

    sum = bit_array[i];
    if (i > 0) {			// special case
	int j;

	j = i - (i & -i);
	i--;
	while (i != j) {		// subtract difference between i and i-1
	    sum -= bit_array[i];
	    i -= i & -i;
	}
    }
    return sum;
}

/**
**	Binary index tree.
*/
static int BitFind(const unsigned *bit_array, int size, unsigned sum)
{
    int idx;
    int bitmask;

    idx = 0;
#if 1
    bitmask = 1 << (sizeof(size) * 8 - 1 - __builtin_clz(size));
#else
    bitmask = size;
    bitmask |= bitmask >> 1;
    bitmask |= bitmask >> 2;
    bitmask |= bitmask >> 4;
    bitmask |= bitmask >> 8;
    bitmask |= bitmask >> 16;
    bitmask = (bitmask >> 1) + 1;
#endif
    while (bitmask) {
	int tmp;

	tmp = idx + bitmask;
	// we need to check, if we allow none 2^n bit arrays
	if (tmp < size) {
	    if (sum == bit_array[tmp]) {	// hit exact
		return tmp;
	    }
	    if (sum > bit_array[tmp]) {	// found smaller
		sum -= bit_array[tmp];
		idx = tmp;
	    }
	}
	bitmask >>= 1;
    }
    return idx;
}

/// @}

// ----------------------------------------------------------------------------

uint8_t Buffer[4096 * 64];
uint64_t EncodeCount;

/**
**	Encode symbol with high probability.
*/
void test1(void)
{
    RangeCoder rc[1];
    int i;
    int n;
    unsigned sum;
    unsigned minsum;
    unsigned maxsum;

    for (n = 0; n < 1000; ++n) {
	unsigned count;

	printf("test1: %d ", n);
	minsum = UINT32_MAX;
	maxsum = 0;
	RangeCoderEncodeInit(rc, n, Buffer, sizeof(Buffer), NULL, NULL);
	for (i = 0; i < 1000; ++i) {
	    RangeCoderEncode(rc, n, 9999, n + 9999);
	}
	count = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += count;

	if (0) {
	    printf("\n");
	    for (i = 0; (unsigned)i < count; ++i) {
		printf("%02x ", Buffer[i]);
	    }
	    printf("\n");
	}

	RangeCoderDecodeInit(rc, &sum, Buffer, sizeof(Buffer), count, NULL,
	    NULL);
	for (i = 0; i < 1000; ++i) {
	    sum = RangeCoderDecode(rc, n + 9999);
	    //printf("%d ", sum);
	    if (sum < minsum) {
		minsum = sum;
	    }
	    if (sum > maxsum) {
		maxsum = sum;
	    }
	    assert(sum > (unsigned)n);
	    assert(sum < (unsigned)n + 9999);
	    RangeCoderDecodeUpdate(rc, n, 9999, n + 9999);
	}
	printf("min %d max %d\r", minsum, maxsum);
    }
}

void test1a(void)
{
    RangeCoder rc[1];
    int i;
    int n;
    unsigned sum;
    unsigned minsum;
    unsigned maxsum;

    for (n = 1; n < 1000; ++n) {
	unsigned count;

	printf("test1a: %d ", n);
	minsum = UINT32_MAX;
	maxsum = 0;
	RangeCoderEncodeInit(rc, n, Buffer, sizeof(Buffer), NULL, NULL);
	for (i = 0; i < 9999; ++i) {
	    RangeCoderEncode(rc, 9999, n, n + 9999);
	}
	count = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += count;

	if (0) {
	    printf("\n");
	    for (i = 0; (unsigned)i < count; ++i) {
		printf("%02x ", Buffer[i]);
	    }
	    printf("\n");
	}

	RangeCoderDecodeInit(rc, &sum, Buffer, sizeof(Buffer), count, NULL,
	    NULL);
	for (i = 0; i < 9999; ++i) {
	    sum = RangeCoderDecode(rc, n + 9999);
	    if (sum < minsum) {
		minsum = sum;
	    }
	    if (sum > maxsum) {
		maxsum = sum;
	    }
	    assert(sum > (unsigned)n);
	    assert(sum < (unsigned)n + 9999);
	    RangeCoderDecodeUpdate(rc, 9999, n, n + 9999);
	}
	printf("min %d max %d\r", minsum, maxsum);
    }
}

void test2()
{
    RangeCoder rc[1];
    volatile unsigned gate_before;
    unsigned frequency[256 + 1];
    volatile unsigned gate_after;
    unsigned sym;
    int i;
    int j;
    int n;
    int k;

    n = 100;
    gate_before = 0x12345678;
    gate_after = 0x87654321;

    for (j = 0; j < 256; j += 16) {
	for (k = 0; (unsigned)k < RangeCoderTotSum - 256;
	    k += (RangeCoderTotSum + (1 << 14) - 1) / (1 << 14)) {
	    unsigned count;

	    printf("test2: %3d/%10d\r", j, k);

	    RangeCoderEncodeInit(rc, 0xAA, Buffer, sizeof(Buffer), NULL, NULL);

	    // all 1 only j + k
	    BitInit(frequency, 256 + 1);
	    for (i = 0; i < 256; ++i) {
		BitUpdate(frequency, 256 + 1, i + 1, 1);
		assert(i == BitFind(frequency, 256 + 1, BitGetSum(frequency,
			    i)));
	    }
	    BitUpdate(frequency, 256 + 1, j + 1, k);
	    assert(BitGetSum(frequency, j) == j);
	    assert(BitGetSum(frequency, j + 1) == j + 1 + k);
	    assert(BitGetSum(frequency, 255) == 255 + k);
	    assert(BitGetSum(frequency, 256) == 256 + k);

	    sym = j;
	    for (i = 0; i < n; ++i) {
		RangeCoderEncode(rc, BitGetSum(frequency, sym),
		    BitGet(frequency, sym + 1), BitGetSum(frequency, 256));
	    }
	    count = RangeCoderEncodeFlush(rc, NULL);
	    EncodeCount += count;

	    if (0) {
		printf("\n");
		for (i = 0; (unsigned)i < count; ++i) {
		    printf("%02x ", Buffer[i]);
		}
		printf("\n");
	    }
	    //if (j == 16) getchar();

	    RangeCoderDecodeInit(rc, &sym, Buffer, sizeof(Buffer), count, NULL,
		NULL);
	    //printf("%02x ", sym);
	    for (i = 0; i < n; ++i) {
		sym = RangeCoderDecode(rc, BitGetSum(frequency, 256));
		//printf("%02x:", sym);
		sym = BitFind(frequency, 256 + 1, sym);

		//printf("=>%02x ", sym);
		//printf("\n");

		assert(sym == (unsigned)j);
		RangeCoderDecodeUpdate(rc, BitGetSum(frequency, sym),
		    BitGet(frequency, sym + 1), BitGetSum(frequency, 256));
	    }
	    //printf("\n");
	}
    }
    assert(gate_before == 0x12345678);
    assert(gate_after == 0x87654321);
}

void test3()
{
    RangeCoder rc[1];
    unsigned frequency[256 + 1];
    unsigned sym;
    int i;
    int n;
    int k;

    printf("			       \r");

    n = 255;
    for (k = 0; (unsigned)k < RangeCoderTotSum - 256;
	k += (RangeCoderTotSum + (1 << 14) - 1) / (1 << 14)) {
	unsigned count;

	printf("test3: %10d\r", k);

	RangeCoderEncodeInit(rc, 0xAA, Buffer, sizeof(Buffer), NULL, NULL);

	BitInit(frequency, 256 + 1);
	for (i = 0; i < 256 + 1; ++i) {
	    BitUpdate(frequency, 256 + 1, i + 1, 1);
	    assert(i == BitFind(frequency, 256 + 1, BitGetSum(frequency, i)));
	}
	BitUpdate(frequency, 256 + 1, 255 + 1, k);

	for (i = 0; i < n; ++i) {
	    sym = i;
	    RangeCoderEncode(rc, BitGetSum(frequency, sym), BitGet(frequency,
		    sym + 1), BitGetSum(frequency, 256));
	}
	count = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += count;

	if (0) {
	    printf("\n");
	    for (i = 0; (unsigned)i < count; ++i) {
		printf("%02x ", Buffer[i]);
	    }
	    printf("\n");
	}

	RangeCoderDecodeInit(rc, &sym, Buffer, sizeof(Buffer), count, NULL,
	    NULL);
	//printf("%02x ", sym);
	for (i = 0; i < n; ++i) {
	    sym = RangeCoderDecode(rc, BitGetSum(frequency, 256));
	    //printf("%02x:", sym);
	    sym = BitFind(frequency, 256 + 1, sym);

	    //printf("=>%02x ", sym);
	    //printf("\n");
	    assert(sym == (unsigned)i);
	    RangeCoderDecodeUpdate(rc, BitGetSum(frequency, sym),
		BitGet(frequency, sym + 1), BitGetSum(frequency, 256));
	}
	//printf("\n");
    }
}

/**
**	Toggle between rarest and most common code.
*/
void test4(void)
{
    RangeCoder rc[1];
    int i;
    int n;
    unsigned sum;

    printf("			       \r");
    for (n = 1; n < 9999; ++n) {
	int count;

	printf("test4: %d\r", n);
	RangeCoderEncodeInit(rc, n, Buffer, sizeof(Buffer), NULL, NULL);
	for (i = 0; i < n; ++i) {
	    RangeCoderEncode(rc, RangeCoderTotSum - 1, 1, RangeCoderTotSum);
	    RangeCoderEncode(rc, 0, RangeCoderTotSum, RangeCoderTotSum);
	}
	count = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += count;

	if (0) {
	    printf("\n");
	    for (i = 0; i < count; ++i) {
		printf("%02x ", Buffer[i]);
	    }
	    printf("\n");
	}

	RangeCoderDecodeInit(rc, &sum, Buffer, sizeof(Buffer), count, NULL,
	    NULL);
	for (i = 0; i < n; ++i) {
	    sum = RangeCoderDecode(rc, RangeCoderTotSum);
	    assert(sum >= RangeCoderTotSum - 1 && sum < RangeCoderTotSum);
	    RangeCoderDecodeUpdate(rc, RangeCoderTotSum - 1, 1,
		RangeCoderTotSum);
	    sum = RangeCoderDecode(rc, RangeCoderTotSum);
	    assert(sum < RangeCoderTotSum);
	    RangeCoderDecodeUpdate(rc, 0, RangeCoderTotSum, RangeCoderTotSum);
	}
    }
}

/**
**	Encode the most uncommon symbol.
**	At all positions.
*/
void test5(void)
{
    RangeCoder rc[1];
    int i;
    int j;
    int n;
    unsigned sum;

    printf("			       \r");
    for (n = 1; n < 100; ++n) {
	for (j = 0; (unsigned)j < RangeCoderTotSum;
	    j += (RangeCoderTotSum + (1 << 14) - 1) / (1 << 14)) {
	    int count;

	    printf("test5: %3d/%10d\r", n, j);

	    RangeCoderEncodeInit(rc, n, Buffer, sizeof(Buffer), NULL, NULL);
	    for (i = 0; i < n; ++i) {
		RangeCoderEncode(rc, j, 1, RangeCoderTotSum);
	    }
	    count = RangeCoderEncodeFlush(rc, NULL);
	    EncodeCount += count;

	    if (0) {
		printf("\n");
		for (i = 0; i < count; ++i) {
		    printf("%02x ", Buffer[i]);
		}
		printf("\n");
	    }

	    RangeCoderDecodeInit(rc, &sum, Buffer, sizeof(Buffer), count, NULL,
		NULL);
	    for (i = 0; i < n; ++i) {
		sum = RangeCoderDecode(rc, RangeCoderTotSum);
		assert(sum == (unsigned)j);
		RangeCoderDecodeUpdate(rc, j, 1, RangeCoderTotSum);
	    }
	}
    }
}

/**
*/
void test6(void)
{
    RangeCoder rc[1];
    int i;
    int n;

    printf("			       \r");
    for (n = 1; n < 100000; ++n) {
	int count;
	unsigned sum;

	RangeCoderEncodeInit(rc, n, NULL, 0, NULL, NULL);

	// 0 1 4686
	for (i = 0; i < 256; ++i) {
	    RangeCoderEncode(rc, 0, 32, 33);
	}
	count = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += count;
#if 0
	for (j = 0; (unsigned)j < RangeCoderTotSum;
	    j += (RangeCoderTotSum + (1 << 14) - 1) / (1 << 14)) {

	    printf("test5: %3d/%10d\r", n, j);

	    for (i = 0; i < n; ++i) {
		RangeCoderEncode(rc, j, 1, RangeCoderTotSum);
	    }

	    if (0) {
		printf("\n");
		for (i = 0; i < count; ++i) {
		    printf("%02x ", Buffer[i]);
		}
		printf("\n");
	    }

	    RangeCoderDecodeInit(rc, &sum, Buffer, sizeof(Buffer), count, NULL,
		NULL);
	    for (i = 0; i < n; ++i) {
		sum = RangeCoderDecode(rc, RangeCoderTotSum);
		assert(sum == (unsigned)j);
		RangeCoderDecodeUpdate(rc, j, 1, RangeCoderTotSum);
	    }
	}
#endif
    }
}

/**
**	Binary search symbol sum.
*/
static inline int findInterval(int *cm, int size, int point)
{
    int index;
    int left;
    int right;
    int cnt;

    left = 0;
    cnt = 0;
    index = -1;
    right = size - 2;
    for (;;) {
	int mid = (right + left) >> 1;

	if (point >= cm[mid] && point < cm[mid + 1]) {
	    index = mid;
	    break;
	}
	if (point >= cm[mid + 1])
	    left = mid + 1;
	else
	    right = mid;
	if (cnt++ >= size)
	    break;
    }
    return index;
}

#include <string.h>

/**
**	Sample file encoder.
*/
void EncodeFile(const char *name_in, const char *name_out)
{
    FILE *file_in;
    FILE *file_out;
    RangeCoder *rc;
    int i;
    int ch;
    unsigned sums[257 + 1];
    unsigned freq[257];
    unsigned count;
    uint8_t *buffer;

    file_in = fopen(name_in, "rb");
    if (!file_in) {
	fprintf(stderr, "can't open input file '%s'\n", name_in);
	return;
    }
    file_out = fopen(name_out, "wb");
    if (!file_out) {
	fprintf(stderr, "can't open output file '%s'\n", name_out);
	fclose(file_in);
	return;
    }

    ch = fgetc(file_in);
    if (ch == EOF) {
	fclose(file_out);
	fclose(file_in);
	return;
    }
    rc = RangeCoderEncodeInit(NULL, ch, Buffer, sizeof(Buffer),
	RangeCoderWriteStdio, file_out);
    if (!rc) {
	fprintf(stderr, "can't initialize range-coder\n");
	fclose(file_out);
	fclose(file_in);
	return;
    }
    // prepare unknown state encode escape as 1
    BitInit(sums, 257 + 1);
    for (i = 0; i < 257; ++i) {
	BitUpdate(sums, 257 + 1, i + 1, 1);
	freq[i] = 1;
    }

    for (;;) {
	// simple order - 1 prediction
	BitUpdate(sums, 257 + 1, ch + 1, 1);
	++freq[ch];
	// rescale if too big
	if (freq[ch] > 127) {
	    BitInit(sums, 257 + 1);
	    for (i = 0; i < 257; ++i) {
		freq[i] = (freq[i] + 1) >> 1;
		BitUpdate(sums, 257 + 1, i + 1, freq[i]);
	    }
	}

	if ((ch = fgetc(file_in)) == EOF) {
	    break;
	}
	assert(ch >= 0);
	RangeCoderEncode(rc, BitGetSum(sums, ch), freq[ch], BitGetSum(sums,
		257));
    }

    // Symbol 256 is our EOF marker
    RangeCoderEncode(rc, BitGetSum(sums, 256), freq[256], BitGetSum(sums,
	    257));
    count = RangeCoderEncodeFlush(rc, &buffer);
    EncodeCount += count;
    // Flush out last buffer
    RangeCoderWriteStdio(file_out, &buffer, &count);

    free(rc);
    fclose(file_out);
    fclose(file_in);
}

/**
**	Sample file decoder.
*/
void DecodeFile(const char *name_in, const char *name_out)
{
    FILE *file_in;
    FILE *file_out;
    RangeCoder *rc;
    int i;
    int ch;
    unsigned sums[257 + 1];
    unsigned freq[257];
    unsigned sym;

    file_in = fopen(name_in, "rb");
    if (!file_in) {
	fprintf(stderr, "can't open input file '%s'\n", name_in);
	return;
    }
    file_out = fopen(name_out, "wb");
    if (!file_out) {
	fprintf(stderr, "can't open output file '%s'\n", name_out);
	fclose(file_in);
	return;
    }

    rc = RangeCoderDecodeInit(NULL, &sym, Buffer, sizeof(Buffer), 0,
	RangeCoderReadStdio, file_in);
    if (!rc) {
	fprintf(stderr, "can't initialize range-coder\n");
	fclose(file_out);
	fclose(file_in);
	return;
    }
    // prepare unknown state encode escape as 1
    BitInit(sums, 257 + 1);
    for (i = 0; i < 257; ++i) {
	BitUpdate(sums, 257 + 1, i + 1, 1);
	freq[i] = 1;
    }

    ch = sym;
    fputc(ch, file_out);
    for (;;) {
	unsigned total;

	// simple order - 1 prediction
	BitUpdate(sums, 257 + 1, ch + 1, 1);
	++freq[ch];
	// rescale if too big
	if (freq[ch] > 127) {
	    BitInit(sums, 257 + 1);
	    for (i = 0; i < 257; ++i) {
		freq[i] = (freq[i] + 1) >> 1;
		BitUpdate(sums, 257 + 1, i + 1, freq[i]);
	    }
	}
	total = BitGetSum(sums, 257);
	sym = RangeCoderDecode(rc, total);
	ch = BitFind(sums, 257 + 1, sym);
	if (ch >= 256) {		// EOF
	    break;
	}
	fputc(ch, file_out);
	RangeCoderDecodeUpdate(rc, BitGetSum(sums, ch), freq[ch], total);
    }

    free(rc);
    fclose(file_out);
    fclose(file_in);
}

#include <time.h>

/**
**	Get ticks in ms.
**
**	@returns ticks in ms,
*/
static inline uint32_t GetMsTicks(void)
{
#ifdef CLOCK_PROCESS_CPUTIME_ID
    struct timespec tspec;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tspec);
    return (tspec.tv_sec * 1000) + (tspec.tv_nsec / (1000 * 1000));
#else
    struct timeval tval;

    if (gettimeofday(&tval, NULL) < 0) {
	return 0;
    }
    return (tval.tv_sec * 1000) + (tval.tv_usec / 1000);
#endif
}

/**
**	Benchmark.
*/
void Bench(void)
{
    RangeCoder rc[1];
    uint32_t start_tick;
    uint32_t end_tick;
    int n;
    int i;

    printf("			       \r");
    printf("bench:\r");
    start_tick = GetMsTicks();
    for (n = 0; n < 10; ++n) {
	RangeCoderEncodeInit(rc, 0xFF, NULL, 0, NULL, NULL);
	for (i = 0; i < 1000; ++i) {
	    RangeCoderEncode(rc, 0, 9999, 9999);
	    RangeCoderEncode(rc, 0, 1, 9999 + 1);
	    RangeCoderEncode(rc, 1, 9999, 9999 + 1);
	    RangeCoderEncode(rc, 1, 1, 9999 + 2);
	}
	i = RangeCoderEncodeFlush(rc, NULL);
	EncodeCount += i;
    }
    end_tick = GetMsTicks();

    printf("\n%d: %.2fs\n", i, (end_tick - start_tick) / 60.0);
}

int main(int argc, const char *argv[])
{
    if (argc > 1) {
	if (argc == 4) {
	    if (!strcmp(argv[1], "-c") || !strcmp(argv[1], "c")) {
		EncodeFile(argv[2], argv[3]);
		return 0;
	    }
	    if (!strcmp(argv[1], "-d") || !strcmp(argv[1], "d")) {
		DecodeFile(argv[2], argv[3]);
		return 0;
	    }
	}
	printf("bad options '-c in out' or '-d in out'\n");
	return 0;
    }

    test1();
    test1a();
    test2();
    test3();
    test4();
    test5();
    //test6();

    Bench();

    // 32: 267793567
    // 64: 545869101
    printf("Encoded %jd bytes\n", (uintmax_t) EncodeCount);
    printf("%d loops max\n", LoopsMax);
    return 0;
}

/// @}
