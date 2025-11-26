#include "openfhe.h"
#include "high-prec-complex.h"

void test_big_fixed_point() {
    BigFixedPoint a = BigFixedPoint(-3.75).scaleTo(128);
    BigFixedPoint b = BigFixedPoint(2.5).scaleTo(128);

    BigFixedPoint c = a + b;
    BigFixedPoint d = a - b;
    BigFixedPoint e = a * b;
    BigFixedPoint f = a / b;

    std::cout << "a: " << a.toString() << " " << a.toBinary() << "\n";
    std::cout << "b: " << b.toString() << " " << b.toBinary() << "\n";
    std::cout << "a + b: " << c.toString() << "\n";
    std::cout << "a - b: " << d.toString() << "\n";
    std::cout << "a * b: " << e.toString() << "\n";
    std::cout << "a / b: " << f.toString() << "\n";
}

void test_zu1() {
    auto b = BigInteger("350551088990192644005587638168552186653");

    auto zu    = getZU();
    auto zuInv = getZUInverse();
    std::cout << zu[0][1].getReal().getValue() << std::endl;
    std::cout << "zu matrix:\n";
    for (size_t i = 0; i != 2; ++i) {
        for (size_t j = 0; j != 4; ++j) {
            std::cout << "zu[" << i << "][" << j << "] = " << zu[i][j].toString() << "\n";
        }
        for (size_t j = zN - 4; j != zN; ++j) {
            std::cout << "zu[" << i << "][" << j << "] = " << zu[i][j].toString() << "\n";
        }
    }
    std::cout << "zuInv matrix:\n";
    for (size_t i = 0; i != 4; ++i) {
        for (size_t j = 0; j != 2; ++j) {
            std::cout << "zuInv[" << i << "][" << j << "] = " << zuInv[i][j].toString() << "\n";
        }
    }
}

void test_zu() {
    std::vector<BigFixedPoint> input = {1, 2, 3, 4, 5, 6, 7, 8};
    input.resize(zN, BigFixedPoint(0, 0, false));
    auto ZU          = getZU();
    auto ZUInv       = getZUInverse();
    auto zuResult    = multU(ZU, input);
    auto zuInvResult = multUInverse(ZUInv, zuResult);

    std::cout << "Input: ";
    for (const auto& val : input) {
        std::cout << std::setprecision(20) << val.toString() << " ";
    }
    std::cout << "\n";

    std::cout << "After ZU and inverse ZU: ";
    for (const auto& val : zuInvResult) {
        std::cout << std::setprecision(20) << val.toString() << " ";
    }
    std::cout << "\n";
}

void test_ru() {
    auto U    = getRU();
    auto UInv = getRUInverse();
    std::cout << "U matrix:\n";
    for (size_t i = 0; i != 4; ++i) {
        for (size_t j = 0; j != 8; ++j) {
            std::cout << "U[" << i << "][" << j << "] = " << U[i][j].toString() << "\n";
        }
    }
    std::cout << "UInv matrix:\n";
    for (size_t i = 0; i != 8; ++i) {
        for (size_t j = 0; j != 4; ++j) {
            std::cout << "UInv[" << i << "][" << j << "] = " << UInv[i][j].toString() << "\n";
        }
    }

    std::vector<BigFixedPoint> input = {1, 2, 3, 4, 5, 6, 7, 8};
    input.resize(rN, BigFixedPoint(0, 0, false));

    auto uResult    = multU(U, input);
    auto uInvResult = multUInverse(UInv, uResult);

    std::cout << "Input: ";
    for (const auto& val : input) {
        std::cout << std::setprecision(20) << val.toString() << " ";
    }
    std::cout << "\n";

    std::cout << "After RU and inverse RU: ";
    for (const auto& val : uInvResult) {
        std::cout << std::setprecision(20) << val.convertToDouble() << " ";
    }
    std::cout << "\n";
}

void test_round() {
    BigFixedPoint a = BigFixedPoint(-3.5).scaleTo(128);
    BigFixedPoint b = BigFixedPoint(2.5).scaleTo(128);
    BigFixedPoint c = BigFixedPoint(3.5).scaleTo(128);
    BigFixedPoint d = BigFixedPoint(-2.5).scaleTo(128);
    BigFixedPoint e = BigFixedPoint(-2.1).scaleTo(128);
    BigFixedPoint f = BigFixedPoint(-0).scaleTo(128);
    BigFixedPoint g = BigFixedPoint(0).scaleTo(128);

    std::cout << "a: " << a.toString() << " rounded: " << a.round().toString() << "\n";
    std::cout << "b: " << b.toString() << " rounded: " << b.round().toString() << "\n";
    std::cout << "c: " << c.toString() << " rounded: " << c.round().toString() << "\n";
    std::cout << "d: " << d.toString() << " rounded: " << d.round().toString() << "\n";
    std::cout << "e: " << e.toString() << " rounded: " << e.round().toString() << "\n";
    std::cout << "f: " << f.toString() << " rounded: " << f.round().toString() << "\n";
    std::cout << "g: " << g.toString() << " rounded: " << g.round().toString() << "\n";
}

void test_ceil() {
    BigFixedPoint a = BigFixedPoint(-3.2).scaleTo(128);
    BigFixedPoint b = BigFixedPoint(2.3).scaleTo(128);
    BigFixedPoint c = BigFixedPoint(3.0).scaleTo(128);
    BigFixedPoint d = BigFixedPoint(-2.7).scaleTo(128);
    BigFixedPoint e = BigFixedPoint(-3.0).scaleTo(128);

    std::cout << "a: " << a.toString() << " ceil: " << a.ceil().toString() << "\n";
    std::cout << "b: " << b.toString() << " ceil: " << b.ceil().toString() << "\n";
    std::cout << "c: " << c.toString() << " ceil: " << c.ceil().toString() << "\n";
    std::cout << "d: " << d.toString() << " ceil: " << d.ceil().toString() << "\n";
    std::cout << "e: " << e.toString() << " ceil: " << e.ceil().toString() << "\n";
}

void test_floor() {
    BigFixedPoint a = BigFixedPoint(-3.8).scaleTo(128);
    BigFixedPoint b = BigFixedPoint(2.7).scaleTo(128);
    BigFixedPoint c = BigFixedPoint(3.0).scaleTo(128);
    BigFixedPoint d = BigFixedPoint(-2.1).scaleTo(128);
    BigFixedPoint e = BigFixedPoint(-3.0).scaleTo(128);

    std::cout << "a: " << a.toString() << " floor: " << a.floor().toString() << "\n";
    std::cout << "b: " << b.toString() << " floor: " << b.floor().toString() << "\n";
    std::cout << "c: " << c.toString() << " floor: " << c.floor().toString() << "\n";
    std::cout << "d: " << d.toString() << " floor: " << d.floor().toString() << "\n";
    std::cout << "e: " << e.toString() << " floor: " << e.floor().toString() << "\n";
}

int main() {
    test_big_fixed_point();
    test_zu1();
    test_zu();
    test_ru();
    test_round();
    test_ceil();
    test_floor();
    return 0;
}