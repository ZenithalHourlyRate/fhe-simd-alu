#include <cassert>
#include "math/z-encode.h"

void test_encodeInZ() {
    auto input   = 255;
    auto encoded = ZPolynomial::encode(input);
    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encoded.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    auto decoded = ZPolynomial::decode(encoded);
    std::cout << "Decoded: " << decoded;
    std::cout << "\n";
}

void test2_encodeInZ() {
    auto input1   = 1;
    auto input2   = 255;
    auto encoded1 = ZPolynomial::encode(input1);
    auto encoded2 = ZPolynomial::encode(input2);
    auto sum      = ZPolynomial::addRaw(encoded1, encoded2);
    std::cout << "Sum 255 + 1 in RE: ";
    for (const auto& val : sum.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << "\n";
    auto rounded = ZPolynomial::roundBigI(sum);
    std::cout << "Sum 255 + 1 in RE Round: ";
    for (const auto& val : rounded.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;
    auto decoded = ZPolynomial::decode(rounded);
    std::cout << "Decoded sum: " << decoded << "\n";
    auto encoded3 = ZPolynomial::encode(input1 + input2);
    std::cout << "Encoded direct 255 + 1 in Z: ";
    for (const auto& val : encoded3.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
}

uint32_t add(uint32_t a, uint32_t b) {
    auto encoded1 = ZPolynomial::encode(a);
    auto encoded2 = ZPolynomial::encode(b);
    auto sum      = ZPolynomial::add(encoded1, encoded2);
    auto decoded  = ZPolynomial::decode(sum);
    return decoded;
}

uint32_t mult(uint32_t a, uint32_t b) {
    auto encoded1 = ZPolynomial::encode(a);
    auto encoded2 = ZPolynomial::encode(b);
    auto res      = ZPolynomial::multiply(encoded1, encoded2);
    auto decoded  = ZPolynomial::decode(res);
    return decoded;
}

void test_add_mult() {
    std::vector<uint32_t> aVec = {1, 2, 3, 33, 255, INT_MAX - 1, INT_MAX, 0xffffffff};
    std::vector<uint32_t> bVec = {1, 2, 3, 255, INT_MAX};

    for (auto a : aVec) {
        for (auto b : bVec) {
            uint32_t addRes  = a + b;
            uint32_t multRes = a * b;
            auto zAddRes     = add(a, b);
            auto zMultRes    = mult(a, b);
            if (addRes != zAddRes) {
                std::cout << "Add Error: a: " << a << " b: " << b << " add: " << addRes << " zAdd: " << zAddRes
                          << std::endl;
            }
            if (multRes != zMultRes) {
                std::cout << "Mult Error: a: " << a << " b: " << b << " mult: " << multRes << " zAdd: " << zMultRes
                          << std::endl;
            }
        }
    }
}

void test3_encodeInZ() {
    uint32_t input1 = 33;
    auto input2     = 57;
    auto encoded1   = ZPolynomial::encode(input1);
    std::cout << "encoded1: " << std::endl;
    for (const auto& val : encoded1.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto encoded2 = ZPolynomial::encode(input2);
    auto mult     = ZPolynomial::multiply(encoded1, encoded2);
    std::cout << "Mult 33 * 57 in Z: ";
    for (const auto& val : mult.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;
    auto decoded = ZPolynomial::decode(mult);
    std::cout << "Decoded 33 * 57 from RE: " << decoded << "\n";
}

void test4_encodeInZ() {
    uint32_t input = -1;
    auto encodedZ  = ZPolynomial::encode(input);

    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encodedZ.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto encodedC = multU(getZU(), encodedZ.getCoefficients());
    //std::cout << "Encoded in C: ";
    //for (const auto& val : encodedC) {
    //    std::cout << val.toHexString() << " ";
    //}
    //std::cout << std::endl;
    auto encodedR = multUInverse(getRUInverse(), encodedC);

    // now back
    auto encCBack = multU(getRU(), encodedR);
    auto encZBack = multUInverse(getZUInverse(), encCBack);
    std::cout << "\nEncoded back in Z: ";
    for (const auto& val : encZBack) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;
    auto decoded = ZPolynomial::decode(encZBack);
    std::cout << "\nDecoded back: " << decoded << "\n";

    // extract error
    auto error = ZPolynomial::extractError(encZBack);
    auto norm  = error.getLog2Norm();
    std::cout << "Error log2 norm: " << norm << std::endl;
    //for (auto e : error.getCoefficients()) {
    //    std::cout << "Error term: " << e.toHexString() << " Norm: " << e.log2Norm() << std::endl;
    //}
    // trunc error
    auto trunc = ZPolynomial::truncError(encZBack);
    std::cout << "Trunc term[0]: " << trunc[0].toHexString() << std::endl;
}

void test_natural() {
    uint32_t input = 255;
    auto encodedZ  = ZPolynomial::encode(input);

    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encodedZ.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto roundNatural = ZPolynomial::roundBigI(encodedZ);
    std::cout << "Rounded Natural in Z: ";
    for (const auto& val : roundNatural.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto multT = ZPolynomial::multiplyRaw(roundNatural, ZPolynomial::getT());
    std::cout << "Rounded Natural * T in Z: ";
    for (const auto& val : multT.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto decoded = ZPolynomial::decode(roundNatural);
    std::cout << "Decoded Rounded Natural * T: " << decoded << "\n";
}

void test_extractI() {
    uint32_t input = 255;
    auto encodedZ  = ZPolynomial::encode(input);
    auto mult      = ZPolynomial::multiplyRaw(ZPolynomial::multiplyRaw(encodedZ, encodedZ), ZPolynomial::getT());

    std::cout << "ExtractI: " << input << "\nEncoded in Z: ";
    for (const auto& val : mult.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto Ipoly = ZPolynomial::extractI(mult);
    std::cout << "Extracted I in Z: ";
    for (const auto& val : Ipoly.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto multNatural = ZPolynomial::roundNOneToOne(mult);
    Ipoly            = ZPolynomial::extractI(multNatural);
    std::cout << "Extracted I in Z after Natural rounding: ";
    for (const auto& val : Ipoly.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;
}

void test_balanced() {
    std::cout << "Test Balanced Encoding\n";
    uint32_t input = 255;
    auto encodedZ  = ZPolynomial::encode(input);

    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encodedZ.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto encodedBalanced = ZPolynomial::toBalancedTInv(encodedZ);
    std::cout << "Encoded Balanced TInv in Z: ";
    for (const auto& val : encodedBalanced.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    auto roundBalanced = ZPolynomial::roundNHalfToHalf(encodedBalanced);
    std::cout << "Rounded Balanced in Z: ";
    for (const auto& val : roundBalanced.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto multT = ZPolynomial::multiplyRaw(roundBalanced, ZPolynomial::getT());
    std::cout << "Rounded Balanced * T in Z: ";
    for (const auto& val : multT.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto roundBalancedToStandard = ZPolynomial::toStandard(roundBalanced);
    std::cout << "Rounded Balanced to Standard in Z: ";
    for (const auto& val : roundBalancedToStandard.getCoefficients()) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto decoded = ZPolynomial::decode(roundBalancedToStandard);
    std::cout << "Decoded Rounded Balanced * T: " << decoded << "\n";
}

int main() {
    test_encodeInZ();
    test2_encodeInZ();
    test3_encodeInZ();
    test_add_mult();
    test4_encodeInZ();
    test_natural();
    test_extractI();
    test_balanced();
    return 0;
}