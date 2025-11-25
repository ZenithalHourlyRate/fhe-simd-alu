#include <cassert>
#include "openfhe.h"
#include "high-prec-complex.h"
#include "z-encode.h"

void test_encodeInZ() {
    auto input   = 255;
    auto encoded = encodeInZ(input);
    std::cout << "Input: " << input << "\nEncoded in Z: ";
    for (const auto& val : encoded) {
        std::cout << val.toBinary(48) << " ";
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
        std::cout << val.toBinary() << " ";
    }
    std::cout << "\n";
    auto rounded = roundInZ(sum);
    std::cout << "Sum 255 + 1 in RE Round: ";
    for (const auto& val : rounded) {
        std::cout << val.toBinary() << " ";
    }
    auto decoded = decodeFromZ(rounded);
    std::cout << "Decoded sum: " << decoded << "\n";
}

void test3_encodeInZ() {
    auto input1      = 33;
    auto input2      = 57;
    auto encoded1    = encodeInZ(input1);
    auto encoded2    = encodeInZ(input2);
    auto multResult  = std::vector<double>(encoded1.size() + encoded2.size() - 1, 0.0);
    auto finalResult = roundInZ(multiplyInRE(encoded1, encoded2));
    auto decoded     = decodeFromZ(finalResult);
    std::cout << "Decoded 33 * 57 from RE: " << decoded << "\n";
}

//void test4_encodeInZ() {
//    auto input    = 3;
//    auto encodedZ = encodeInZ(input);
//    auto encodedC = multiplyByZUComplex(encodedZ);
//    auto encodedR = multiplyByUInverseComplex(encodedC);
//    std::cout << "Input: " << input << "\nEncoded in RE: ";
//    for (const auto& val : encodedR) {
//        std::cout << std::setprecision(20) << val << " ";
//    }
//
//    // now back
//    auto encCBack = multiplyByUComplex(encodedR);
//    auto encZBack = multiplyByZUInverse(encCBack);
//    std::cout << "\nEncoded back in RE: ";
//    for (const auto& val : encZBack) {
//        std::cout << std::setprecision(20) << val << " ";
//    }
//    auto decoded = decodeFromRE(encZBack);
//    std::cout << "\nDecoded back: " << decoded << "\n";
//}

int main() {
    test_encodeInZ();
    test2_encodeInZ();
    test3_encodeInZ();
    //test4_encodeInZ();
    return 0;
}