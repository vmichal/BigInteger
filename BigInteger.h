
// BigInteger (non)standard header
// Vojtìch Michal, 21.3.2018; 

//28.3.2018 pøidán kód umožòující benchmarkování metod
//31.3.2018 dokonèen a debugnut	algoritmus na odmocòování
//25.8.2018 pøechod na C++17
//konec srpna provedeny dílèí optimalizace

#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <vector> 
#include <string>
#include <cstdint>
#include <type_traits>


#define UNIT_BIT_COUNT sizeof(unit)*CHAR_BIT

class BigInteger {

public:
	using unit = std::uint32_t;
	struct msbFirst_t {};
	struct lsbFirst_t {};
	static constexpr inline msbFirst_t msbFirst = { };
	static constexpr inline lsbFirst_t lsbFirst = { };

	template<typename T>
	struct is_order_tag : std::conditional_t<
		std::is_same_v<T, lsbFirst_t> || std::is_same_v<T, msbFirst_t>
		, std::true_type, std::false_type> {};
	template<typename T>
	static constexpr inline bool is_order_tag_v = is_order_tag<T>::value;

private:
	std::vector<unit> cislice_;
	int8_t sign_ = 0;

	std::vector<uint8_t> convert(int8_t base) const;
	inline void checkInternalState();

public:
	/////////////////////////////
	// Konstruktory, assignment
	/////////////////////////////
	BigInteger(unit cislo);
	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	BigInteger(T cislo);

	template<typename Cont, typename = std::enable_if_t<std::is_same_v<
		std::remove_reference_t<Cont>, std::vector<unit>>>>
		BigInteger(Cont&&, int8_t sign = 1);														

	template<typename Iter, typename Order,
		typename = typename std::iterator_traits<Iter>::value_type,
		typename = std::enable_if_t<is_order_tag_v<Order>>>
		BigInteger(Iter begin, Iter end, Order o);

	BigInteger(const std::string& str, int8_t radix = 10);

	BigInteger(const BigInteger& rhs) = default;
	BigInteger& operator=(const BigInteger& rhs) = default;
	BigInteger(BigInteger&& rhs) = default;
	BigInteger& operator=(BigInteger&& rhs) = default;
	BigInteger() = default;


	/////////////////////////////
	// Pøístup k jednotlivým èíslùm
	/////////////////////////////
	unit& operator[](size_t index);
	const unit& operator[](size_t index) const;
	unit& at(size_t index);
	const unit& at(size_t index) const;


	/////////////////////////////
	// Size, bitCount
	/////////////////////////////
	///<summary>
	///Vrátí poèet 32bitových integerù, které tvoøí èíslo (délka vectoru èísel).
	///</summary>
	size_t size() const { return cislice_.size(); }
	size_t bitCount() const;
	size_t highBitsCount() const;
	explicit operator bool() const {
		return !isZero();
	}

	/////////////////////////////
	// Bitové operace
	/////////////////////////////
	void setBit(std::size_t index);
	void clearBit(std::size_t index);
	void toggleBit(std::size_t index);
	bool isBitSet(std::size_t index) const;

	/////////////////////////////
	// Bitová logika
	/////////////////////////////
	BigInteger operator&(const BigInteger& rhs) const;
	BigInteger& operator&=(const BigInteger& rhs);

	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	BigInteger operator&(T rhs) const;
	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	BigInteger& operator&=(T rhs);
	
	BigInteger operator|(const BigInteger& rhs) const;
	BigInteger& operator|=(const BigInteger& rhs);
	
	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	BigInteger operator|(T rhs) const;
	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	BigInteger& operator|=(T rhs);



	/////////////////////////////
	// Bit shift
	/////////////////////////////
	BigInteger operator<<(std::size_t pocet) const;
	BigInteger operator>>(std::size_t pocet) const;
	BigInteger& operator<<=(std::size_t pocet);
	BigInteger& operator>>=(std::size_t pocet);

	/////////////////////////////
	// Aritmetické operátory
	/////////////////////////////
	BigInteger operator+(const BigInteger &druhy) const;
	BigInteger operator-(const BigInteger &druhy) const;
	BigInteger operator*(const BigInteger &druhy) const;
	BigInteger operator/(const BigInteger &druhy) const;
	BigInteger operator%(const BigInteger &druhy) const;

	/////////////////////////////
	// Aritmetické operátory a assignment
	/////////////////////////////
	BigInteger& operator+=(const BigInteger &druhy);
	BigInteger& operator-=(const BigInteger &druhy);
	BigInteger& operator*=(const BigInteger &druhy);
	BigInteger& operator/=(const BigInteger &druhy);
	BigInteger& operator%=(const BigInteger &druhy);


	/////////////////////////////
	// Inkrementace, dekrementace
	/////////////////////////////
	BigInteger operator++(int);
	BigInteger operator--(int);
	BigInteger& operator++();
	BigInteger& operator--();

	/////////////////////////////
	// Mocnìní
	/////////////////////////////
	BigInteger pow(unit mocnina) const;
	BigInteger sqrt() const;

	/////////////////////////////
	// Is metody (kladne, zero, liche)
	/////////////////////////////
	bool isZero() const { return !sign_; }
	bool isKladne() const { return sign_ == 1; }
	bool isZaporne() const { return sign_ == -1; }
	bool isSude() const { return !isLiche(); }
	bool isLiche() const { return cislice_[0] & 0x1; }


	/////////////////////////////
	// Comparison, equality
	/////////////////////////////
	bool operator<(const BigInteger& druhy) const;
	bool operator>(const BigInteger& druhy) const;
	bool operator==(const BigInteger& druhy) const;
	bool operator!=(const BigInteger& druhy) const;
	bool operator>=(const BigInteger& druhy) const;
	bool operator<=(const BigInteger& druhy) const;


	/////////////////////////////
	// Aritmetické metody, absolutní hodnoty
	/////////////////////////////
	static BigInteger add(const BigInteger& a, const BigInteger& b);
	static BigInteger sub(const BigInteger& a, const BigInteger& b);
	static BigInteger multiply(const BigInteger& a, const BigInteger& b);
	static BigInteger divide(const BigInteger& delenec_param, const BigInteger& delitel_param, BigInteger* zbytek = nullptr);
	static BigInteger zbytek(const BigInteger& delenec_param, const BigInteger& delitel_param);
	static BigInteger abs(const BigInteger& cislo);
	static BigInteger zapornaAbs(const BigInteger& cislo);
	static BigInteger opacnaHodnota(const BigInteger& cislo);

	/////////////////////////////
	// Convert funkce, print
	/////////////////////////////
	std::string printDec() const;
	std::string printHex() const;
	std::string printBin() const;
	std::string print(int8_t base) const;
	std::string printExp(unit presnost = 0) const;
	int64_t valueOf() const;

	/////////////////////////////
	// Input/Output
	/////////////////////////////
	friend std::ostream& operator<<(std::ostream& stream, const BigInteger& a);
	friend std::istream& operator>>(std::istream& stream, BigInteger& a);
};
#endif //BIGINTEGER_H ^^  