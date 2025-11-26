#include <cassert>
#include "openfhe.h"
#include "high-prec-complex.h"
#include "z-encode.h"

void test_encodeInZ() {
    auto input   = 255;
    auto encoded = encodeInZ(input);
    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encoded) {
        std::cout << val.toHexString() << " ";
    }
    auto decoded = decodeFromZ(encoded);
    std::cout << "Decoded: " << decoded;
    std::cout << "\n";
}

void test2_encodeInZ() {
    auto input1   = 1;
    auto input2   = 255;
    auto encoded1 = encodeInZ(input1);
    auto encoded2 = encodeInZ(input2);
    std::vector<BigFixedPoint> sum(encoded1.size(), 0);
    for (size_t i = 0; i != encoded1.size(); ++i) {
        sum[i] = encoded1[i] + encoded2[i];
    }
    std::cout << "Sum 255 + 1 in RE: ";
    for (const auto& val : sum) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << "\n";
    auto rounded = roundInZ(sum);
    std::cout << "Sum 255 + 1 in RE Round: ";
    for (const auto& val : rounded) {
        std::cout << val.toHexString() << " ";
    }
    auto decoded = decodeFromZ(rounded);
    std::cout << "Decoded sum: " << decoded << "\n";
    auto encoded3 = encodeInZ(input1 + input2);
    std::cout << "Encoded direct 255 + 1 in Z: ";
    for (const auto& val : encoded3) {
        std::cout << val.toHexString() << " ";
    }
}

uint32_t add(uint32_t a, uint32_t b) {
    auto encoded1 = encodeInZ(a);
    auto encoded2 = encodeInZ(b);
    std::vector<BigFixedPoint> sum(encoded1.size(), 0);
    for (size_t i = 0; i != encoded1.size(); ++i) {
        sum[i] = encoded1[i] + encoded2[i];
    }
    //std::cout << "Sum 255 + 1 in RE: ";
    //for (const auto& val : sum) {
    //    std::cout << val.toHexString() << " ";
    //}
    //std::cout << "\n";
    auto rounded = roundInZ(sum);
    //std::cout << "Sum 255 + 1 in RE Round: ";
    //for (const auto& val : rounded) {
    //    std::cout << val.toHexString() << " ";
    //}
    auto decoded = decodeFromZ(rounded);
    //std::cout << "Decoded sum: " << decoded << "\n";
    //auto encoded3 = encodeInZ(input1 + input2);
    //std::cout << "Encoded direct 255 + 1 in Z: ";
    //for (const auto& val : encoded3) {
    //    std::cout << val.toHexString() << " ";
    //}
    return decoded;
}

uint32_t mult(uint32_t a, uint32_t b) {
    auto encoded1 = encodeInZ(a);
    //std::cout << "encoded1: " << std::endl;
    //for (const auto& val : encoded1) {
    //    std::cout << val.toHexString() << " ";
    //}
    //std::cout << std::endl;

    auto encoded2 = encodeInZ(b);
    auto mult     = multiplyInRE(encoded1, encoded2);
    //std::cout << "Mult 33 * 57 in RE: ";
    //for (const auto& val : mult) {
    //    std::cout << val.toHexString() << " ";
    //}
    //std::cout << std::endl;
    auto finalResult = roundInZ(mult);
    auto decoded     = decodeFromZ(finalResult);
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
    auto encoded1   = encodeInZ(input1);
    std::cout << "encoded1: " << std::endl;
    for (const auto& val : encoded1) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto encoded2 = encodeInZ(input2);
    auto mult     = multiplyInRE(encoded1, encoded2);
    std::cout << "Mult 33 * 57 in RE: ";
    for (const auto& val : mult) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;
    auto finalResult = roundInZ(mult);
    auto decoded     = decodeFromZ(finalResult);
    std::cout << "Decoded 33 * 57 from RE: " << decoded << "\n";
}

void test4_encodeInZ() {
    auto input    = 7;
    auto encodedZ = encodeInZ(input);

    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encodedZ) {
        std::cout << val.toHexString() << " ";
    }
    std::cout << std::endl;

    auto encodedC = multU(getZU(), encodedZ);
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
    auto decoded = decodeFromZ(encZBack);
    std::cout << "\nDecoded back: " << decoded << "\n";
}

int main() {
    //test_encodeInZ();
    //test2_encodeInZ();
    //test3_encodeInZ();
    test_add_mult();
    test4_encodeInZ();
    return 0;
}