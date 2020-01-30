#include "BigInteger.h"	  									 
#include <exception>
#include <cctype>
#include <algorithm>
#include <functional>


#define BIGINTEGER_DEBUG  
#ifdef BIGINTEGER_DEBUG
#define SCOPE_TIMER Timer timer(__FUNCTION__, __LINE__); 
#else
#define SCOPE_TIMER
#endif

#ifdef BIGINTEGER_DEBUG	 
#include <iostream>
#include <chrono>
#endif

namespace {

#ifdef BIGINTEGER_DEBUG

	class Timer {
		static unsigned next_id_;
		std::chrono::time_point<std::chrono::steady_clock> start_;
		unsigned id_, line_;
		const char *functionName_;
	public:
		Timer(const char* functionName, unsigned line) {
			id_ = next_id_++;
			functionName_ = functionName;
			std::cout << "Benchmark id " << id_ << " in " << functionName << " started\n";
			start_ = std::chrono::high_resolution_clock::now();
		}

		~Timer() {
			std::chrono::milliseconds dur = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_);
			std::cout << "Benchmark " << id_ << " of function " << functionName_ << " on line " << line_ << " ended in " << dur.count() << "ms\n";
		}
	};
	unsigned Timer::next_id_ = 0;
#endif

	namespace implementation {
		template<typename Iter>
		std::vector<unit> genFromBool(Iter begin, Iter end) {
			static_assert(std::is_same_v<std::iterator_traits<Iter>::value_type, bool>, "Sem se mìl dostat jen bool, blbèe");
			SCOPE_TIMER;
			std::vector<unit> vec;
			while (begin != end) {
				unit tmp = 0;
				for (unit bit = 1; begin != end && bit; bit <<= 1, ++begin)
					if (*begin)
						tmp |= bit;
				vec.push_back(tmp);
			}
			return vec;
		}

		template<typename Iter>
		std::vector<unit> genFromGreater(Iter begin, Iter end) {
			static_assert(sizeof(std::iterator_traits<Iter>::value_type) > sizeof(BigInteger::unit));
			SCOPE_TIMER;
			std::vector<unit> vec;
			constexpr int multiple = sizeof(std::iterator_traits<Iter>::value_type) / sizeof(BigInteger::unit);
			for (; begin != end; ++begin) {
				std::iterator_traits<Iter>::value_type value = *begin;
				for (int i = 0; i < multiple; ++i) {
					vec.push_back(static_cast<unit>(value));
					value >>= UNIT_BIT_COUNT;
				}
			}
			return vec;
		}

		template<typename Iter>
		std::vector<unit> genFromSmaller(Iter begin, Iter end) {
			static_assert(sizeof(std::iterator_traits<Iter>::value_type) < sizeof(BigInteger::unit));
			SCOPE_TIMER;
			std::vector<unit> vec;
			constexpr int multiple = sizeof(unit) / sizeof(std::iterator_traits<Iter>::value_type);
			for (;;) {
				unit tmp = 0;
				for (int i = 0; i < multiple && begin != end; ++i, ++begin)
					tmp |= *begin << sizeof(*begin)*CHAR_BIT;
				vec.push_back(tmp);
			}
			return vec;
		}
	}

	template<typename Iter, typename = std::iterator_traits<Iter>::value_type>
	std::vector<unit> generateVector(Iter begin, Iter end, BigInteger::lsbFirst_t) {
		SCOPE_TIMER;
		using Val = std::iterator_traits<Iter>::value_type;
		if constexpr (sizeof(Val) == sizeof(unit))
			return std::vector<unit>(begin, end);
		if constexpr (std::is_same_v<Val, bool>)
			return implementation::genFromBool(begin, end);
		if constexpr (sizeof(Val) > sizeof(unit))
			return implementation::genFromGreater(begin, end);
		if constexpr (sizeof(Val) < sizeof(unit))
			return implementation::genFromSmaller(begin, end);
	}

	template<typename Iter, typename = std::iterator_traits<Iter>::value_type>
	std::vector<unit> generateVector(Iter begin, Iter end, BigInteger::msbFirst_t) {
		SCOPE_TIMER;
		std::vector<std::iterator_traits<Iter>::value_type> temporary;
		if constexpr (std::is_same_v<std::iterator_traits<Iter>::iterator_category, std::random_access_iterator_tag>)
			temporary.reserve(static_cast<std::size_t>(std::distance(begin, end));
		std::copy(begin, end, std::back_inserter(temporary));
		return generateVector(temporary.rbegin(), temporary.rend(), BigInteger::lsbFirst);
	}
}

BigInteger::BigInteger(unit cislo)
	: sign_(cislo > 0) {
	if (sign_)
		cislice_.push_back(cislo);
}

template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
BigInteger::BigInteger(T cislo)
	: sign_(cislo > 0 ? 1 : cislo < 0 ? -1 : 0) {
	SCOPE_TIMER;
	uint64_t c = cislo < 0 ? -cislo : cislo;
	while (c) {
		cislice_.push_back(static_cast<unit>(c));
		c >>= 32;
	}
	checkInternalState();
}

template<typename Cont, typename = std::enable_if_t<std::is_same_v<std::remove_reference_t<Cont>, std::vector<BigInteger::unit>>>>
BigInteger::BigInteger(Cont&& vec, int8_t sign)
	: sign_(sign), cislice_(vec) {
	SCOPE_TIMER;
	checkInternalState();
}


template<typename Iter, typename Order,
	typename = typename std::iterator_traits<Iter>::value_type,
	typename = std::enable_if_t<BigInteger::is_order_tag_v<Order>>>
	BigInteger::BigInteger(Iter begin, Iter end, Order o)
	: sign_(1),
	cislice_(generateVector(begin, end, o)) {
	checkInternalState();
}

BigInteger::BigInteger(const std::string& str, int8_t radix) {
	SCOPE_TIMER;
	if (radix < 2)
		throw std::invalid_argument("Radix below 2");
	auto iterator = str.begin(), end = str.end();
	while (iterator != end && *iterator == '0')
		++iterator;
	if (end == iterator)
		return;
	if (*iterator == '-') {
		sign_ = -1;
		++iterator;
	}
	else
		sign_ = 1;
	const unsigned pocetZnakuNa32bitu = static_cast<unsigned>(std::log10(static_cast<uint64_t>(UINT32_MAX) + 1) / std::log10(radix)) + 1;
	cislice_.reserve(std::distance(iterator, end) / pocetZnakuNa32bitu);	   //rezervuj dost místa pro pøepoèet všech znakù
	for (; iterator != end; ++iterator) {
		if (*iterator == '\'' || *iterator == '_')
			continue;
		unsigned char c = static_cast<unsigned char>(*iterator);
		uint64_t carry = isdigit(c) ? c - '0' : islower(c) ? c - 'a' + 10 : c - 'A' + 10;
		if (carry >= radix) //neplatná èíslice v zápisu èísla
			break;
		for (auto &i : cislice_) {
			carry += static_cast<uint64_t>(i) * radix;
			i = static_cast<unit>(carry);
			carry >>= 0x20;
		}
		if (carry)
			cislice_.push_back(static_cast<unit>(carry));
	}
	if (cislice_.size())
		checkInternalState();
	else
		sign_ = 0;
}

BigInteger::unit& BigInteger::operator[](std::size_t index) {
	SCOPE_TIMER;
	return cislice_[index];
}

const BigInteger::unit& BigInteger::operator[](std::size_t index) const {
	SCOPE_TIMER;
	return cislice_[index];
}

BigInteger::unit& BigInteger::at(std::size_t index) {
	SCOPE_TIMER;
	return cislice_.at(index);
}

const BigInteger::unit& BigInteger::at(std::size_t index) const {
	SCOPE_TIMER;
	return cislice_.at(index);
}

///<summary>
///Odstarní z vectoru všechny pøebyteèné nuly a v pøípadì potøeby upraví znaménko
///</summary>
void BigInteger::checkInternalState() {
	SCOPE_TIMER;
	while (cislice_.size() && !cislice_.back())
		cislice_.pop_back();
	if (cislice_.empty())
		sign_ = 0;
	else if (!sign_)
		throw std::exception("vnitøní stav objektu byl narušen");
}

std::size_t BigInteger::bitCount() const {
	SCOPE_TIMER;
	if (isZero())
		return 0;
	std::size_t count = (size() - 1) << 5;
	for (unit highestByte = cislice_.back(); highestByte; count++, highestByte >>= 1);
	return count;
}

std::size_t BigInteger::highBitsCount() const {
	SCOPE_TIMER;
	std::size_t count = 0;
	for (auto i : cislice_)
		for (; i; i ^= i & -i, ++count);
	/*for (auto iter = cislice_.begin(); iter != cislice_.end(); ++iter)
		for (int i = 0; i < 32; i++)
			if (*iter & 1 << i)			LUL, já byl tak hloupej
				++count; */
	return count;
}

void BigInteger::setBit(std::size_t index) {
	SCOPE_TIMER;
	if (isZero())
		sign_ = 1;
	std::size_t indexCisla = index >> 5;
	if (indexCisla >= size())
		cislice_.resize(indexCisla + 1);
	cislice_[indexCisla] |= 1 << (index & 0x1f);
}

void BigInteger::clearBit(std::size_t index) {
	SCOPE_TIMER;
	if (index < bitCount()) {
		cislice_[index >> 5] &= ~(1 << (index & 0x1f));
		checkInternalState();
	}
}

void BigInteger::toggleBit(std::size_t index) {
	SCOPE_TIMER;
	if (index >= bitCount())
		return setBit(index);
	cislice_[index >> 5] ^= 1 << (index & 0x1f);
	checkInternalState();
}

bool BigInteger::isBitSet(std::size_t index) const {
	SCOPE_TIMER;
	if (index >= bitCount())
		return false;
	return cislice_[index >> 5] & 1 << (index & 0x1f);
}

BigInteger BigInteger::operator&(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) &= rhs;
}

BigInteger& BigInteger::operator&=(const BigInteger& rhs) {
	SCOPE_TIMER;
	if (isZero())
		return *this;
	if (rhs.isZero())
		return *this = BigInteger();
	if (std::size_t rhsSize = rhs.cislice_.size(); rhs < cislice_.size())
		cislice_.resize(rhsSize);
	std::transform(cislice_.begin(), cislice_.end(), rhs.cislice_.begin(), cislice_.begin(), std::bit_and<unit>());
	return *this;
}

template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
BigInteger BigInteger::operator&(T rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) &= rhs;
}

template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
BigInteger& BigInteger::operator&=(T rhs) {
	SCOPE_TIMER;
	if (isZero())
		return *this;
	if (!rhs)
		return *this = BigInteger();
	std::vector<unit> novy;
	for (auto iter = cislice_.begin(), end = cislice_.end(); rhs && iter != end; ++iter, rhs >>= UNIT_BIT_COUNT)
		novy.push_back(*iter & static_cast<unit>(rhs));
	cislice_.swap(novy);
	return *this;
}

BigInteger BigInteger::operator|(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) |= rhs;
}

BigInteger& BigInteger::operator|=(const BigInteger& rhs) {
	SCOPE_TIMER;
	if (rhs.isZero())
		return *this;
	if (isZero())
		return *this = rhs;
	if (cislice_.capacity() < rhs.size()) {
		std::vector<unit> novy;
		novy.reserve(rhs.size());
		std::transform(cislice_.begin(), cislice_.end(), rhs.cislice_.begin(), std::back_inserter(novy), std::bit_or<unit>());
		std::copy(rhs.cislice_.begin() + novy.size(), rhs.cislice_.end(), std::back_inserter(novy));
		cislice_.swap(novy);
	}
	else
		std::transform(rhs.cislice_.begin(), rhs.cislice_.end(), cislice_.begin(), cislice_.begin(), std::bit_or<unit>());
	return *this;
}

template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
BigInteger BigInteger::operator|(T rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) |= rhs;
}

template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
BigInteger& BigInteger::operator|=(T rhs) {
	SCOPE_TIMER;
	if (isZero())
		return *this = rhs;
	if (cislice_.size() * sizeof(unit) < sizeof(T))
		cislice_.resize(sizeof(T) / sizeof(unit));
	for (auto iter = cislice_.begin(); rhs; rhs >>= UNIT_BIT_COUNT, ++iter)
		*iter |= static_cast<unit>(rhs);
	return *this;
}


BigInteger BigInteger::operator<<(std::size_t pocet) const {
	SCOPE_TIMER;
	return BigInteger(*this) <<= pocet;
}

BigInteger& BigInteger::operator<<=(std::size_t pocet) {
	SCOPE_TIMER;
	if (isZero())
		return *this;
	std::vector<unit> novy;
	novy.reserve(cislice_.size() + (pocet >> 5) + 1);
	novy.resize(pocet >> 5);
	pocet &= 0x1f;
	uint64_t carry = 0;
	std::transform(cislice_.begin(), cislice_.end(), std::back_inserter(novy), [pocet, &carry](unit u) mutable -> unit {
		carry = carry >> UNIT_BIT_COUNT | static_cast<uint64_t>(u) << pocet;
		return static_cast<unit>(carry);
	});
	for (carry >>= UNIT_BIT_COUNT; carry; carry >>= UNIT_BIT_COUNT)
		novy.push_back(static_cast<unit>(carry));
	cislice_.swap(novy);
	return *this;
}

BigInteger BigInteger::operator>>(std::size_t pocet) const {
	SCOPE_TIMER;
	return BigInteger(*this) >>= pocet;
}

BigInteger& BigInteger::operator>>=(std::size_t pocet) {
	SCOPE_TIMER;
	if (isZero())
		return *this;
	if (pocet >= bitCount())
		return *this = BigInteger();
	auto iterator = cislice_.begin() + (pocet >> 5);
	pocet &= 0x1f;
	uint64_t carry = *iterator >> pocet;
	auto newEnd = std::transform(++iterator, cislice_.end(), cislice_.begin(), [&carry, pocet](unit u) mutable -> unit {
		carry = carry >> sizeof(unit)*CHAR_BIT | static_cast<uint64_t>(u) << pocet;
		return static_cast<unit>(carry);
	});
	cislice_.resize(std::distance(cislice_.begin(), newEnd));
	for (carry >>= UNIT_BIT_COUNT, carry; carry >>= UNIT_BIT_COUNT)
		cislice_.push_back(static_cast<unit>(carry));
	/*

	auto iterator = cislice_.begin() + (pocet >> 5 < cislice_.size() ? pocet >> 5 : cislice_.size()); //posune iterator o (pocet/32) bitù nebo na konec vectoru
	if (iterator == cislice_.end()) //bylo-li posunutí vìtší, poté vracíme nulu
		return BigInteger();
	pocet &= 0x1f; //totéž co modulo 32
	if (!pocet)
		return BigInteger(std::vector<unit>(iterator, cislice_.end()), sign_);
	std::vector<unit> novy;
	uint64_t carry = *iterator++ >> pocet;
	while (iterator != cislice_.end()) {
		carry |= uint64_t(*iterator++) << (0x20 - pocet);
		novy.push_back(carry & 0xFFFF'FFFF);
		carry >>= 32;
	}
	if (carry)
		novy.push_back(carry & 0xFFFF'FFFF);
	return BigInteger(std::move(novy), sign_);*/
}

BigInteger BigInteger::operator+(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) += rhs;
}

BigInteger BigInteger::operator-(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) -= rhs;
}

BigInteger BigInteger::operator*(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) *= rhs;
}

BigInteger BigInteger::operator/(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) /= rhs;
}

BigInteger BigInteger::operator%(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return BigInteger(*this) %= rhs;
}

BigInteger BigInteger::operator++(int) {
	SCOPE_TIMER;
	BigInteger res(*this);
	++*this;
	return res;
}

BigInteger BigInteger::operator--(int) {
	SCOPE_TIMER;
	BigInteger res(*this);
	--*this;
	return res;
}

BigInteger& BigInteger::operator++() {
	SCOPE_TIMER;
	return *this += 1;
}

BigInteger& BigInteger::operator--() {
	SCOPE_TIMER;
	return *this -= 1;
}

BigInteger& BigInteger::operator+=(const BigInteger& rhs) {
	SCOPE_TIMER;
	if (isZero())
		return *this = rhs;
	if (rhs.isZero())
		return *this;
	if (std::size_t rhsSize = rhs.size(); cislice_.capacity() <= rhsSize) {
		std::vector<unit> novy;
		novy.reserve(rhs.size() + 1);
		std::copy(rhs.cislice_.begin(), rhs.cislice_.end(), std::back_inserter(novy));
		uint64_t carry = 0;
		std::transform(cislice_.begin(), cislice_.end(), novy.begin(), novy.begin(), [carry](unit a, unit b) mutable -> unit {
			carry = (carry >> UNIT_BIT_COUNT) + a + b;
			return static_cast<unit>(carry);
		});
		cislice_.swap(novy);
	}
	else {
		cislice_.resize(rhsSize);
		uint64_t carry = 0;
		std::transform(cislice_.begin(), cislice_.end(), rhs.cislice_.begin(), cislice_.begin(), [&carry](unit a, unit b) mutable ->unit {
			carry = (carry >> UNIT_BIT_COUNT) + a + b;
			return static_cast<unit>(carry);
		});
		for (carry >>= UNIT_BIT_COUNT; carry; carry >>= UNIT_BIT_COUNT)
			cislice_.push_back(static_cast<unit>(carry));
	}
	return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& rhs) {
	SCOPE_TIMER;

}

BigInteger& BigInteger::operator*=(const BigInteger& rhs) {
	SCOPE_TIMER;

}

BigInteger& BigInteger::operator/=(const BigInteger& rhs) {
	SCOPE_TIMER;

}

BigInteger& BigInteger::operator%=(const BigInteger& rhs) {
	SCOPE_TIMER;
}


BigInteger BigInteger::pow(unit mocnina) const {
	SCOPE_TIMER;
	switch (mocnina) {
		case 0:
			return 1;
		case 1:
			return *this;
		case 2:
			return *this**this;
	}
	mocnina -= 2;
	BigInteger vysledek = *this**this;

	while (mocnina--)
		vysledek *= *this;
	return vysledek;
}

///<summary>
/// Vrátí nejvìtší celé èíslo menší nebo rovno odmocninì
///<summary>
BigInteger BigInteger::sqrt() const {
	SCOPE_TIMER;
	if (isZero())
		return 0;
	/*if (bitCount() < 64) {
	return BigInteger(std::sqrt(valueOf()));
	}*/
	int shift = (bitCount() & 0x1f) - 2;
	if (shift & 1)
		shift++;
	auto iterator = cislice_.rbegin();
	BigInteger zbytek = *iterator >> shift & 0x3, odmocnina = std::sqrtl(zbytek.valueOf());
	zbytek -= odmocnina.pow(2);
	for (; iterator != cislice_.rend(); iterator++) {
		while ((shift -= 2) >= 0) {
			zbytek = (zbytek << 2) + (*iterator >> shift & 0x3);
			unit cislice = zbytek >= (odmocnina << 2 | 1);
			if (cislice)
				zbytek -= (odmocnina << 2) | 1;
			odmocnina = odmocnina << 1 | cislice;
		}
		shift = 32;
	}
	return odmocnina;
}

bool BigInteger::operator<(const BigInteger& rhs) const {
	SCOPE_TIMER;
	if (sign_ != rhs.sign_)
		return sign_ < rhs.sign_;
	if (bitCount() != rhs.bitCount())
		return isKladne() == (bitCount() < rhs.bitCount()); //vrátí true pro (sign_ = 1 && this<rhs) nebo (sign_ = -1 && this>rhs
	auto[thisIter, otherIter] = std::mismatch(cislice_.rbegin(), cislice_.rend(), rhs.cislice_.rbegin());
	if (thisIter == cislice_.rend())
		return false;
	return isKladne() == (*thisIter < *otherIter);
}

bool BigInteger::operator>(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return rhs < *this;
}

bool BigInteger::operator==(const BigInteger & rhs) const {
	SCOPE_TIMER;
	return sign_ == rhs.sign_
		&& bitCount() == rhs.bitCount()
		&& cislice_.end() == std::mismatch(cislice_.begin(), cislice_.end(), rhs.cislice_.begin()).first;
}

bool BigInteger::operator!=(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return !(*this == rhs);
}

bool BigInteger::operator<=(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return !(*this > rhs);
}

bool BigInteger::operator>=(const BigInteger& rhs) const {
	SCOPE_TIMER;
	return !(*this < rhs);
}

BigInteger BigInteger::add(const BigInteger& lhs, const BigInteger& rhs) {
	SCOPE_TIMER;
	if (lhs.isZero())
		return rhs;
	if (rhs.isZero())
		return lhs;
	if (lhs.isKladne() && rhs.isKladne()) {
		uint64_t carry = 0;
		const std::vector<unit>& kratsi = lhs.size() < rhs.size() ? lhs.cislice_ : rhs.cislice_;
		std::vector<unit> vysledek(lhs.size() < rhs.size() ? rhs.cislice_ : lhs.cislice_);
		auto k = kratsi.begin();
		auto v = vysledek.begin();

		while (k != kratsi.end()) {
			carry += uint64_t(*k++) + *v;
			*v++ = carry & 0xFFFF'FFFF;
			carry >>= 32;
		}
		while (v != vysledek.end() && carry) {
			carry += *v;
			*v++ = carry & 0xFFFF'FFFF;
			carry >>= 32;
		}
		if (carry)
			vysledek.push_back(carry);
		return BigInteger(std::move(vysledek), 1);
	}
	if (lhs.isZaporne() && rhs.isZaporne())
		return zapornaAbs(abs(lhs) + abs(rhs));

	if (lhs.isZaporne() && rhs.isKladne())
		return rhs - abs(lhs);

	if (lhs.isKladne() && rhs.isZaporne())
		return lhs - abs(rhs);
}

BigInteger BigInteger::sub(const BigInteger& mensenec, const BigInteger& mensitel) {
	SCOPE_TIMER;
	if (mensitel.isZero())
		return mensenec;
	if (mensenec.isZero())
		return opacnaHodnota(mensitel);
	if (mensenec.isKladne() && mensitel.isKladne()) {
		if (mensenec >= mensitel) {
			std::vector<unit> vysledek(mensenec.cislice_);
			auto delsi = vysledek.begin();
			auto kratsi = mensitel.cislice_.begin();
			int64_t carry = 0;
			while (kratsi != mensitel.cislice_.end()) {
				carry += int64_t(*delsi) - *kratsi++;
				if (carry >= 0) {
					*delsi++ = carry;
					carry = 0;
				}
				else {
					*delsi++ = carry + 0x1'0000'0000;
					carry = -1;
				}
			}
			while (carry) {
				carry += int64_t(*delsi);
				if (carry >= 0) {
					*delsi++ = carry;
					carry = 0;
				}
				else {
					*delsi++ = carry + 0xFFFF'FFFF;
					carry = -1;
				}
			}
			return BigInteger(std::move(vysledek), 1);
		}
		return zapornaAbs(mensitel - mensenec);
	}
	if (mensenec.isKladne() && mensitel.isZaporne())
		return mensenec + abs(mensitel);

	if (mensenec.isZaporne() && mensitel.isKladne())
		return zapornaAbs(abs(mensenec) + mensitel);

	if (mensenec.isZaporne() && mensitel.isZaporne())
		return opacnaHodnota(abs(mensenec) - abs(mensitel));
}

BigInteger BigInteger::multiply(const BigInteger& a, const BigInteger& b) {
	SCOPE_TIMER;
	if (a.isZero() || b.isZero())
		return BigInteger();
	std::vector<unit> novy(a.cislice_.size() + b.cislice_.size() + 1);
	const std::vector<unit>& delsi = a.size() < b.size() ? b.cislice_ : a.cislice_,
		&kratsi = a.size() < b.size() ? a.cislice_ : b.cislice_;
	uint64_t carry, offset = 0;
	for (std::size_t i = 0; i < kratsi.size(); i++) {
		carry = 0;
		for (std::size_t j = 0; j < delsi.size(); j++) {
			carry += uint64_t(kratsi[i]) * delsi[j] + novy[i + j];
			novy[i + j] = carry & 0xffff'ffff;
			carry >>= 32;
		}
		if (carry)
			novy[i + delsi.size()] = carry & 0xffff'ffff;
	}
	return BigInteger(std::move(novy), a.sign_ == b.sign_ ? 1 : -1);
}

BigInteger BigInteger::divide(const BigInteger& delenec_param, const BigInteger& delitel_param, BigInteger* zbytek) {
	SCOPE_TIMER;
	if (delitel_param.isZero())
		throw std::exception("Není možné dìlit nulou");
	if (delitel_param == 1)
		return delenec_param;
	if (delitel_param == -1)
		return opacnaHodnota(delenec_param);
	bool trebaMazatZbytek = false;
	if (!zbytek) {
		trebaMazatZbytek = true;
		zbytek = new BigInteger(abs(delenec_param));
	}
	else
		new(zbytek) BigInteger(abs(delenec_param));
	BigInteger podil;
	BigInteger* delitel = new BigInteger(1);
	unit mocnina = 0;
	std::vector<BigInteger> mocniny;
	mocniny.reserve(delitel_param.bitCount() / delitel_param.bitCount());
	while (*zbytek >= *delitel) {
		mocniny.push_back(*delitel);
		mocnina++;
		*delitel *= abs(delitel_param);
	}
	delete delitel;
	if (mocniny.empty()) {
		if (trebaMazatZbytek)
			delete zbytek;
		else
			zbytek->sign_ = delenec_param.sign_ == delitel_param.sign_ ? delenec_param.sign_ : -delenec_param.sign_;
		return BigInteger();
	}
	delitel = &mocniny.back();
	while (--mocnina) {
		unit miraZanoreni = 0;
		while (*zbytek >= *delitel) {
			miraZanoreni++;
			*delitel <<= 1;
		}
		while (miraZanoreni--) {
			*delitel >>= 1;
			if (*zbytek >= *delitel) {
				*zbytek -= *delitel;
				podil += mocniny[mocnina - 1] << miraZanoreni;
			}
		}
		delitel--;
	}
	podil.sign_ = zbytek->sign_ = delitel_param.sign_ == delenec_param.sign_ ? delenec_param.sign_ : -delenec_param.sign_;
	if (trebaMazatZbytek)
		delete zbytek;
	return podil;
}

BigInteger BigInteger::zbytek(const BigInteger& delenec, const BigInteger& delitel) {
	SCOPE_TIMER;
	BigInteger zbytek;
	divide(delenec, delitel, &zbytek);
	return zbytek;
}

BigInteger BigInteger::abs(const BigInteger& cislo) {
	SCOPE_TIMER;
	if (cislo.isZero())
		return BigInteger();
	BigInteger novy(cislo);
	novy.sign_ = 1;
	return novy;
}

BigInteger BigInteger::zapornaAbs(const BigInteger& cislo) {
	SCOPE_TIMER;
	if (cislo.isZero())
		return BigInteger();
	BigInteger res = cislo;
	res.sign_ = -1;
	return res;
}

BigInteger BigInteger::opacnaHodnota(const BigInteger& cislo) {
	SCOPE_TIMER;
	if (cislo.isZero())
		return BigInteger();
	BigInteger res(cislo);
	res.sign_ = cislo.isKladne() ? -1 : 1;
	return res;
}

std::string BigInteger::printDec() const {
	SCOPE_TIMER;
	return print(10);
}

std::string BigInteger::printHex() const {
	SCOPE_TIMER;
	if (isZero())
		return "0"
	if (sign_ == -1)
		str.push_back('-');
	str.push_back('0');
	str.push_back('x');
	unit prvni = *cislice_.rbegin();
	bool jsouCisla = false;
	for (int i = 7; i >= 0; i--) {
		int nibble = prvni >> (4 * i) & 0xf;
		if (nibble || jsouCisla) {
			str.push_back(nibble < 10 ? nibble + '0' : nibble + 'A' - 10);
			jsouCisla = true;
		}
		if (jsouCisla && i % 2 == 0)
			str.push_back('_');
	}
	for (auto iterator = cislice_.rbegin() + 1; iterator != cislice_.rend(); iterator++) {
		unit c = *iterator;
		for (int i = 7; i >= 0; i--) {
			int nibble = c >> (4 * i) & 0xf;
			str.push_back(nibble < 10 ? nibble + '0' : nibble + 'A' - 10);
			if (i % 2 == 0)
				str.push_back('_');
		}
	}
	str.pop_back();
	return str;
}

std::string BigInteger::printBin() const {
	SCOPE_TIMER;
	std::string str;
	if (!sign_) {
		str.push_back('0');
		return str;
	}
	if (sign_ == -1)
		str.push_back('-');
	str.push_back('0');
	str.push_back('b');
	unit prvni = *cislice_.rbegin();
	bool jsouCisla = false;
	for (int i = 31; i >= 0; i--) {
		int8_t bit = prvni >> i & 1;
		if (jsouCisla)
			str.push_back(bit + '0');
		else if (bit) {
			str.push_back(bit + '0');
			jsouCisla = true;
		}
		if (jsouCisla && i % 4 == 0)
			str.push_back('_');
	}
	for (auto iterator = cislice_.rbegin() + 1; iterator != cislice_.rend(); iterator++) {
		unit c = *iterator;
		for (int i = 31; i >= 0; i--) {
			int8_t bit = c >> i & 1;
			str.push_back(bit + '0');
			if (i % 4 == 0)
				str.push_back('_');
		}
	}
	str.pop_back();
	return str;
}

std::string BigInteger::print(int8_t base) const {
	SCOPE_TIMER;
	if (base == 2)
		return printBin();
	if (base == 16)
		return printHex();
	std::vector<uint8_t> cisla(convert(base));
	std::string str;
	if (sign_ == -1)
		str.push_back('-');
	for (std::vector<uint8_t>::const_iterator i = cisla.begin(); i != cisla.end(); i++) {
		str.push_back(*i + '0');
		if ((cisla.end() - i - 1) % 3 == 0)
			str.push_back('.');
	}
	str.pop_back();
	return str;
}

std::string BigInteger::printExp(unit presnost) const {
	SCOPE_TIMER;
	if (sign_ == 0)
		return std::string("0");
	std::vector<uint8_t> cisla(convert(10));
	std::string str;
	std::vector<uint8_t>::const_iterator iterator = cisla.begin();
	if (sign_ == -1)
		str.push_back('-');
	str.push_back(*iterator++ + '0');
	str.push_back(',');
	int8_t rad = 0;
	if (presnost) {
		while (iterator != cisla.end() && presnost--) {
			str.push_back(*iterator++ + '0');
			if (++rad % 3 == 0) {
				str.push_back('.');
			}
		}
		if (iterator != cisla.end())
			if (*iterator > 4)
				str.back()++;
	}
	else
		while (iterator != cisla.end()) {
			str.push_back(*iterator++ + '0');
			if (++rad % 3 == 0) {
				str.push_back('.');
			}
		}
	str.append(" x 10^");
	str.append(std::to_string(cisla.size() - 1));
	return str;
}

///<summary>
///Vrátí posledních 63 bitù èísla jako int64_t hodnotu
///</summary>
int64_t BigInteger::valueOf() const {
	SCOPE_TIMER;
	if (!sign_)
		return 0;
	if (cislice_.size() == 1)
		return cislice_.at(0) * sign_;
	int64_t result = (uint64_t(cislice_.at(1)) << 32 | cislice_.at(0)) & 0x7fff'ffff'ffff'ffff;
	if (sign_ == -1)
		result = ~result + 1;
	return result;
}

///<summary>
///Vrátí èíslo jako vector bajtù reprezentujících èíslice jeho zápisu v soustavì o daném základu
///<summary>
std::vector<uint8_t> BigInteger::convert(int8_t base) const {
	SCOPE_TIMER;
	if (!sign_)
		return std::vector<uint8_t> { 0 };
	std::vector<uint8_t> novy;
	for (auto iterator = cislice_.rbegin(); iterator != cislice_.rend(); iterator++) {
		uint64_t carry = *iterator;
		for (uint8_t & i : novy) {
			carry += uint64_t(i) << 32;
			i = carry % base;
			carry /= base;
		}
		while (carry) {
			novy.push_back(carry % base);
			carry /= base;
		}
	}
	std::reverse(novy.begin(), novy.end());
	return novy;
}

std::ostream & operator<<(std::ostream & stream, const BigInteger& a) {
	SCOPE_TIMER;
	if (stream.flags() & std::ostream::hex)
		return stream << a.printHex();
	if (stream.flags() & std::ostream::binary)
		return stream << a.printBin();
	return stream << a.printDec();
}

std::istream& operator>>(std::istream& stream, BigInteger& a) {
	SCOPE_TIMER;
	std::string str;
	char znak;
	while (isalpha(znak = stream.get()) || isdigit(znak)) {
		str.push_back(znak);
	}
	new(&a) BigInteger(std::move(str));
	return stream;
}