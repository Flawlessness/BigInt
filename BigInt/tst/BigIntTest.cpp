#include "gtest/gtest.h"
#include "../src/BigInt.h"

TEST(BigIntTest, default_constructor_test)
{
    BigInt a;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, int_constructor_max_test)
{
    BigInt a(1234567890);
    EXPECT_EQ(std :: string (a), "1234567890");
}

TEST(BigIntTest, int_constructor_negative_test)
{
    BigInt a(-1234567890);
    EXPECT_EQ(std :: string (a), "-1234567890");
}

TEST(BigIntTest, int_constructor_positive_test)
{
    BigInt a(+1234567890);
    EXPECT_EQ(std :: string (a), "1234567890");
}

TEST(BigIntTest, int_constructor_null_test)
{
    BigInt a(-0);
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, string_constructor_negative_test)
{
    BigInt a("-1234567890");
    EXPECT_EQ(std :: string (a), "-1234567890");
}

TEST(BigIntTest, string_constructor_positive_test)
{
    BigInt a("+1234567890");
    EXPECT_EQ(std :: string (a), "1234567890");
}

TEST(BigIntTest, string_constructor_big_test)
{
    BigInt a("472347895734089573489571894239428357348751897590123482387423545235");
    EXPECT_EQ(std :: string (a), "472347895734089573489571894239428357348751897590123482387423545235");
}

TEST(BigIntTest, string_constructor_small_test)
{
    BigInt a("-546");
    EXPECT_EQ(std :: string (a), "-546");
}

TEST(BigIntTest, string_constructor_null_test_1)
{
    BigInt a("-0");
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, string_constructor_null_test_2)
{
    BigInt a("-0123");
    EXPECT_EQ(std :: string (a), "-123");
}

TEST(BigIntTest, string_constructor_null_test_3)
{
    BigInt a("-0000000000123");
    EXPECT_EQ(std :: string (a), "-123");
}

TEST(BigIntTest, string_constructor_null_test_4)
{
    BigInt a("-000000123");
    EXPECT_EQ(std :: string (a), "-123");
}

TEST(BigIntTest, string_constructor_exception_test_1)
{
    EXPECT_THROW(BigInt a("2132g"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_2)
{
    EXPECT_THROW(BigInt a("21+32"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_3)
{
    EXPECT_THROW(BigInt a(""), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_4)
{
    EXPECT_THROW(BigInt a("-+1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_5)
{
    EXPECT_THROW(BigInt a("+-1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_6)
{
    EXPECT_THROW(BigInt a("+-----1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_7)
{
    EXPECT_THROW(BigInt a("-+++++1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_8)
{
    EXPECT_THROW(BigInt a("-------1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_9)
{
    EXPECT_THROW(BigInt a("+++++++1"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_10)
{
    EXPECT_THROW(BigInt a("-+"), std :: invalid_argument);
}

TEST(BigIntTest, string_constructor_exception_test_11)
{
    EXPECT_THROW(BigInt a("+-"), std :: invalid_argument);
}

TEST(BigIntTest, bigint_constructor_positive_test)
{
    BigInt a("+1947395729547");
    BigInt b(a);
    EXPECT_EQ(std :: string (a), std :: string (b));
}

TEST(BigIntTest, bigint_constructor_negative_test)
{
    BigInt a("-1947395729547");
    BigInt b(a);
    EXPECT_EQ(std :: string (a), std :: string (b));
}

TEST(BigIntTest, operator_assignment_test_1)
{
    BigInt a("-57238573546546345");
    a = a;
    EXPECT_EQ(std :: string (a), "-57238573546546345");
}

TEST(BigIntTest, operator_assignment_test_2)
{
    BigInt a("-57238573546546345"), b;
    b = a;
    EXPECT_EQ(std :: string (b), "-57238573546546345");
}

TEST(BigIntTest, operator_not_negative_small_test)
{
    const int expected_result = -173684;
    BigInt a("-173684");
    EXPECT_EQ(static_cast<int>(~a), ~expected_result);
}

TEST(BigIntTest, operator_not_negative_big_test)
{
    BigInt a("-93475936485663483434563485623865238465346775676456346456546455633423295793475802374523495874");
    a = ~a;
    EXPECT_EQ(std :: string (a), "93475936485663483434563485623865238465346775676456346456546455633423295793475802374523495873");
}

TEST(BigIntTest, operator_not_positive_small_test)
{
    const int expected_result = 173684;
    BigInt a("173684");
    EXPECT_EQ(static_cast<int>(~a), ~expected_result);
}

TEST(BigIntTest, operator_not_positive_big_test)
{
    BigInt a("826946154463785622378462349234785693478562347856384756347856347865782652234234");
    a = ~a;
    EXPECT_EQ(std :: string (a), "-826946154463785622378462349234785693478562347856384756347856347865782652234235");
}

TEST(BigIntTest, operator_not_null_test_1)
{
    const int expected_result = 0;
    BigInt a("0");
    EXPECT_EQ(static_cast<int>(~a), ~expected_result);
}

TEST(BigIntTest, operator_not_null_test_2)
{
    const int expected_result = -1;
    BigInt a("-1");
    EXPECT_EQ(static_cast<int>(~a), ~expected_result);
}

TEST(BigIntTest, operator_prefix_increment_test)
{
    BigInt a("623492347856934785"), b;
    b = ++a;
    EXPECT_EQ(std :: string (b), "623492347856934786");
    EXPECT_EQ(std :: string (a), "623492347856934786");
}

TEST(BigIntTest, operator_postfix_increment_test)
{
    BigInt a("623492347856934785"), b;
    b = a++;
    EXPECT_EQ(std :: string (b), "623492347856934785");
    EXPECT_EQ(std :: string (a), "623492347856934786");
}

TEST(BigIntTest, operator_prefix_decrement_test)
{
    BigInt a("623492347856934785"), b;
    b = --a;
    EXPECT_EQ(std :: string (b), "623492347856934784");
    EXPECT_EQ(std :: string (a), "623492347856934784");
}

TEST(BigIntTest, operator_postfix_decrement_test)
{
    BigInt a("623492347856934785"), b;
    b = a--;
    EXPECT_EQ(std :: string (b), "623492347856934785");
    EXPECT_EQ(std :: string (a), "623492347856934784");
}

TEST(BigIntTest, operator_addition_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a += b;
    EXPECT_EQ(std :: string (a), "623492347861259019");
}

TEST(BigIntTest, operator_addition_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a += b;
    EXPECT_EQ(std :: string (a), "859038539875935201839");
}

TEST(BigIntTest, operator_addition_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a += b;
    EXPECT_EQ(std :: string (a), "859038539875933954855");
}

TEST(BigIntTest, operator_addition_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a += b;
    EXPECT_EQ(std :: string (a), "-859038539875933954855");
}

TEST(BigIntTest, operator_addition_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a += b;
    EXPECT_EQ(std :: string (a), "859038539875933954855");
}

TEST(BigIntTest, operator_addition_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a += b;
    EXPECT_EQ(std :: string (a), "-859038539875933954855");
}

TEST(BigIntTest, operator_addition_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a += b;
    EXPECT_EQ(std :: string (a), "-859038539875935201839");
}

TEST(BigIntTest, operator_addition_assign_test_8)
{
    BigInt a("-645687346798774"), b("-547457566798879");
    a += b;
    EXPECT_EQ(std :: string (a), "-1193144913597653");
}

TEST(BigIntTest, operator_addition_assign_test_9)
{
    BigInt a("2222222222"), b("2222222222");
    a += b;
    EXPECT_EQ(std :: string (a), "4444444444");
}

TEST(BigIntTest, operator_addition_assign_test_10)
{
    BigInt a("-348729375897238947298123463758978347589273947587346326363"), b("6723461278346985734256347895617268365785645");
    a += b;
    EXPECT_EQ(std :: string (a), "-348729375897232223836845116773244091241378330318980540718");
}

TEST(BigIntTest, operator_addition_assign_test_11)
{
    BigInt a("-123"), b("123");
    a += b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_addition_assign_test_12)
{
    BigInt a("123"), b("-123");
    a += b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_addition_assign_test_13)
{
    BigInt a("999999999"), b("999999900");
    a += b;
    EXPECT_EQ(std :: string (a), "1999999899");
}

TEST(BigIntTest, operator_addition_assign_test_14)
{
    BigInt a("46578345634657834657834652345234534"), b("-46578345634657834657834652345234534");
    a += b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_addition_assign_test_15)
{
    BigInt a("-623492"), b("-859038539875934578347");
    a += b;
    EXPECT_EQ(std :: string (a), "-859038539875935201839");
}

TEST(BigIntTest, operator_subtraction_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a -= b;
    EXPECT_EQ(std :: string (a), "623492347852610551");
}

TEST(BigIntTest, operator_subtraction_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a -= b;
    EXPECT_EQ(std :: string (a), "-859038539875933954855");
}

TEST(BigIntTest, operator_subtraction_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a -= b;
    EXPECT_EQ(std :: string (a), "-859038539875935201839");
}

TEST(BigIntTest, operator_subtraction_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a -= b;
    EXPECT_EQ(std :: string (a), "859038539875935201839");
}

TEST(BigIntTest, operator_subtraction_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a -= b;
    EXPECT_EQ(std :: string (a), "859038539875935201839");
}

TEST(BigIntTest, operator_subtraction_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a -= b;
    EXPECT_EQ(std :: string (a), "-859038539875935201839");
}

TEST(BigIntTest, operator_subtraction_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a -= b;
    EXPECT_EQ(std :: string (a), "-859038539875933954855");
}

TEST(BigIntTest, operator_subtraction_assign_test_8)
{
    BigInt a("-645687346798774"), b("-547457566798879");
    a -= b;
    EXPECT_EQ(std :: string (a), "-98229779999895");
}

TEST(BigIntTest, operator_subtraction_assign_test_9)
{
    BigInt a("2222222222"), b("2222222222");
    a -= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_subtraction_assign_test_10)
{
    BigInt a("-348729375897238947298123463758978347589273947587346326363"), b("6723461278346985734256347895617268365785645");
    a -= b;
    EXPECT_EQ(std :: string (a), "-348729375897245670759401810744712603937169564855712112008");
}

TEST(BigIntTest, operator_subtraction_assign_test_11)
{
    BigInt a("123"), b("123");
    a -= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_subtraction_assign_test_12)
{
    BigInt a("-123"), b("-123");
    a -= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_subtraction_assign_test_13)
{
    BigInt a("-9999999999"), b("-9991111111111");
    a -= b;
    EXPECT_EQ(std :: string (a), "9981111111112");
}

TEST(BigIntTest, operator_subtraction_assign_test_14)
{
    BigInt a("-989898989898998"), b("-878732424124235451234342424344234");
    a -= b;
    EXPECT_EQ(std :: string (a), "878732424124235450244443434445236");
}

TEST(BigIntTest, operator_subtraction_assign_test_15)
{
    BigInt a("-623492"), b("-859038539875934578347");
    a -= b;
    EXPECT_EQ(std :: string (a), "859038539875933954855");
}

TEST(BigIntTest, operator_multipl_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a *= b;
    EXPECT_EQ(std :: string (a), "2696126809342784533079690");
}

TEST(BigIntTest, operator_multipl_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a *= b;
    EXPECT_EQ(std :: string (a), "535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a *= b;
    EXPECT_EQ(std :: string (a), "-535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a *= b;
    EXPECT_EQ(std :: string (a), "-535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a *= b;
    EXPECT_EQ(std :: string (a), "-535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a *= b;
    EXPECT_EQ(std :: string (a), "-535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a *= b;
    EXPECT_EQ(std :: string (a), "535603657304326202122727724");
}

TEST(BigIntTest, operator_multipl_assign_test_8)
{
    BigInt a("-645687346798774"), b("0");
    a *= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_multipl_assign_test_9)
{
    BigInt a("0"), b("2222222");
    a *= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_multipl_assign_test_10)
{
    BigInt a("-0"), b("0");
    a *= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_multipl_assign_test_11)
{
    BigInt a("-348729375897238947298123463758978347589273947587346326363"), b("6723461278346985734256347895617268365785645");
    a *= b;
    EXPECT_EQ(std :: string (a), "-2344668455467196687824877236999534484165366006690638916855618194493381029897250148000106111070459135");
}

TEST(BigIntTest, operator_multipl_assign_test_12)
{
    BigInt a("-69729834753906732984534768732759874353434753024573894753894753897458934758934759873689767546974568974"), b("-49783547590346589723467853489753489753894523475634753489738967348957348957234956345623489673895734895634");
    a *= b;
    EXPECT_EQ(std :: string (a), "3471398546938119424101825047499489503130828761790928232193533714811552193068507387340437271312431368483399165595782542703760576714657411026724300132759995181340283437595817207864981135098191234378024459516");
}

TEST(BigIntTest, operator_division_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a /= b;
    EXPECT_EQ(std :: string (a), "144185617118");
}

TEST(BigIntTest, operator_division_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a /= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_division_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a /= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_division_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a /= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_division_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a /= b;
    EXPECT_EQ(std :: string (a), "-1377785985828101");
}

TEST(BigIntTest, operator_division_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a /= b;
    EXPECT_EQ(std :: string (a), "-1377785985828101");
}

TEST(BigIntTest, operator_division_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a /= b;
    EXPECT_EQ(std :: string (a), "1377785985828101");
}

TEST(BigIntTest, operator_division_assign_test_8)
{
    BigInt a("0"), b("-859038539875934578347");
    a /= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_division_assign_test_9)
{
    BigInt a("39875934578347"), b("0");
    EXPECT_THROW(a /= b, std :: invalid_argument);
}

TEST(BigIntTest, operator_division_assign_test_10)
{
    BigInt a("-59734569726975697734534563465397453475460694668976948689732567345897349583945379045239862379"), b("45");
    a /= b;
    EXPECT_EQ(std :: string (a), "-1327434882821682171878545854786610077232459881532821081994057052131052212976563978783108052");
}

TEST(BigIntTest, operator_division_assign_test_11)
{
    BigInt a("-123879127423805734897567384658793533333333333333460534753453475628374534"), b("-573495749578");
    a /= b;
    EXPECT_EQ(std :: string (a), "216007054132416346139351656450113069460900187396263030431657");
}

TEST(BigIntTest, operator_xor_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a ^= b;
    EXPECT_EQ(std :: string (a), "623492347861019659");
}

TEST(BigIntTest, operator_xor_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a ^= b;
    EXPECT_EQ(std :: string (a), "859038539875935135023");
}

TEST(BigIntTest, operator_xor_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a ^= b;
    EXPECT_EQ(std :: string (a), "-859038539875935135017");
}

TEST(BigIntTest, operator_xor_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a ^= b;
    EXPECT_EQ(std :: string (a), "-859038539875935135023");
}

TEST(BigIntTest, operator_xor_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a ^= b;
    EXPECT_EQ(std :: string (a), "-859038539875935135017");
}

TEST(BigIntTest, operator_xor_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a ^= b;
    EXPECT_EQ(std :: string (a), "-859038539875935135023");
}

TEST(BigIntTest, operator_xor_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a ^= b;
    EXPECT_EQ(std :: string (a), "859038539875935135017");
}

TEST(BigIntTest, operator_xor_assign_test_8)
{
    BigInt a("0"), b("-859038539875934578347");
    a ^= b;
    EXPECT_EQ(std :: string (a), "-859038539875934578347");
}

TEST(BigIntTest, operator_xor_assign_test_9)
{
    BigInt a("39875934578347"), b("0");
    a ^= b;
    EXPECT_EQ(std :: string (a), "39875934578347");
}

TEST(BigIntTest, operator_xor_assign_test_10)
{
    BigInt a("39875934578347"), b("-0");
    a ^= b;
    EXPECT_EQ(std :: string (a), "39875934578347");
}

TEST(BigIntTest, operator_xor_assign_test_11)
{
    BigInt a("8590385398759345783471234123454356723846364873462837432423534"), b("43756348795638475634785634785634785634786537846578346534534");
    a ^= b;
    EXPECT_EQ(std :: string (a), "8627864627191194941290207350413384561214080664166099856863976");
}

TEST(BigIntTest, operator_xor_assign_test_12)
{
    BigInt a("4123454356723846364873"), b("56347856347856347856347");
    a ^= b;
    EXPECT_EQ(std :: string (a), "52851647976750009008914");
}

TEST(BigIntTest, operator_mod_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a %= b;
    EXPECT_EQ(std :: string (a), "4297173");
}

TEST(BigIntTest, operator_mod_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a %= b;
    EXPECT_EQ(std :: string (a), "623492");
}

TEST(BigIntTest, operator_mod_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a %= b;
    EXPECT_EQ(std :: string (a), "859038539875933954855");
}

TEST(BigIntTest, operator_mod_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a %= b;
    EXPECT_EQ(std :: string (a), "623492");
}

TEST(BigIntTest, operator_mod_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a %= b;
    EXPECT_EQ(std :: string (a), "229655");
}

TEST(BigIntTest, operator_mod_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a %= b;
    EXPECT_EQ(std :: string (a), "393837");
}

TEST(BigIntTest, operator_mod_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a %= b;
    EXPECT_EQ(std :: string (a), "393837");
}

TEST(BigIntTest, operator_mod_assign_test_8)
{
    BigInt a("0"), b("-859038539875934578347");
    a %= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_mod_assign_test_9)
{
    BigInt a("39875934578347"), b("0");
    EXPECT_THROW(a %= b, std::invalid_argument);
}

TEST(BigIntTest, operator_mod_assign_test_10)
{
    BigInt a("8590385398759345783471234123454356723846364873462837432423534"), b("786");
    a %= b;
    EXPECT_EQ(std :: string (a), "758");
}


TEST(BigIntTest, operator_and_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a &= b;
    EXPECT_EQ(std :: string (a), "119680");
}

TEST(BigIntTest, operator_and_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a &= b;
    EXPECT_EQ(std :: string (a), "33408");
}

TEST(BigIntTest, operator_and_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a &= b;
    EXPECT_EQ(std :: string (a), "859038539875934544936");
}

TEST(BigIntTest, operator_and_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a &= b;
    EXPECT_EQ(std :: string (a), "590084");
}

TEST(BigIntTest, operator_and_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a &= b;
    EXPECT_EQ(std :: string (a), "859038539875934544936");
}

TEST(BigIntTest, operator_and_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a &= b;
    EXPECT_EQ(std :: string (a), "590084");
}

TEST(BigIntTest, operator_and_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a &= b;
    EXPECT_EQ(std :: string (a), "-859038539875935168428");
}

TEST(BigIntTest, operator_and_assign_test_8)
{
    BigInt a("0"), b("-859038539875934578347");
    a &= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_and_assign_test_9)
{
    BigInt a("6485729576"), b("-6485729576");
    a &= b;
    EXPECT_EQ(std :: string (a), "8");
}

TEST(BigIntTest, operator_and_assign_test_10)
{
    BigInt a("8972348503453849573489753578093745893745893745897358455463464564324234534"), b("-5346898274923569345634785683745637846573455345565465127982534543");
    a &= b;
    EXPECT_EQ(std :: string (a), "8972348498272679583808603437692235795808751169499800551363825801202167840");
}

TEST(BigIntTest, operator_or_assign_test_1)
{
    BigInt a("623492347856934785"), b("4324234");
    a |= b;
    EXPECT_EQ(std :: string (a), "623492347861139339");
}

TEST(BigIntTest, operator_or_assign_test_2)
{
    BigInt a("623492"), b("859038539875934578347");
    a |= b;
    EXPECT_EQ(std :: string (a), "859038539875935168431");
}

TEST(BigIntTest, operator_or_assign_test_3)
{
    BigInt a("-623492"), b("859038539875934578347");
    a |= b;
    EXPECT_EQ(std :: string (a), "-590081");
}

TEST(BigIntTest, operator_or_assign_test_4)
{
    BigInt a("623492"), b("-859038539875934578347");
    a |= b;
    EXPECT_EQ(std :: string (a), "-859038539875934544939");
}

TEST(BigIntTest, operator_or_assign_test_5)
{
    BigInt a("859038539875934578347"), b("-623492");
    a |= b;
    EXPECT_EQ(std :: string (a), "-590081");
}

TEST(BigIntTest, operator_or_assign_test_6)
{
    BigInt a("-859038539875934578347"), b("623492");
    a |= b;
    EXPECT_EQ(std :: string (a), "-859038539875934544939");
}

TEST(BigIntTest, operator_or_assign_test_7)
{
    BigInt a("-859038539875934578347"), b("-623492");
    a |= b;
    EXPECT_EQ(std :: string (a), "-33411");
}

TEST(BigIntTest, operator_or_assign_test_8)
{
    BigInt a("0"), b("-859038539875934578347");
    a |= b;
    EXPECT_EQ(std :: string (a), "-859038539875934578347");
}

TEST(BigIntTest, operator_or_assign_test_9)
{
    BigInt a("6485729576"), b("-6485729576");
    a |= b;
    EXPECT_EQ(std :: string (a), "-8");
}

TEST(BigIntTest, operator_or_assign_test_10)
{
    BigInt a("8972348503453849573489753578093745893745893745897358455463464564324234534"), b("-5346898274923569345634785683745637846573455345565465127982534543");
    a |= b;
    EXPECT_EQ(std :: string (a), "-165728285242419205233275585808495270175897441465826364860467849");
}

TEST(BigIntTest, operator_unary_plus_test_1)
{
    BigInt a("546665634534");
    a = +a;
    EXPECT_EQ(std :: string (a), "546665634534");
}

TEST(BigIntTest, operator_unary_plus_test_2)
{
    BigInt a("-546665634534");
    a = +a;
    EXPECT_EQ(std :: string (a), "-546665634534");
}

TEST(BigIntTest, operator_unary_plus_test_3)
{
    BigInt a("-0");
    a = +a;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_unary_minus_test_1)
{
    BigInt a("546665634534");
    a = -a;
    EXPECT_EQ(std :: string (a), "-546665634534");
}

TEST(BigIntTest, operator_unary_minus_test_2)
{
    BigInt a("-546665634534");
    a = -a;
    EXPECT_EQ(std :: string (a), "546665634534");
}

TEST(BigIntTest, operator_unary_minus_test_3)
{
    BigInt a("0");
    a = -a;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_equality_test_1)
{
    BigInt a("0"), b("0");
    a = a == b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_equality_test_2)
{
    BigInt a("1"), b("0");
    a = a == b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_equality_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a == b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_equality_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a == b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_equality_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a == b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_equality_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a == b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_equality_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a == b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_inequality_test_1)
{
    BigInt a("0"), b("0");
    a = a != b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_inequality_test_2)
{
    BigInt a("1"), b("0");
    a = a != b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_inequality_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a != b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_inequality_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a != b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_inequality_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a != b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_inequality_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a != b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_inequality_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a != b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_test_1)
{
    BigInt a("0"), b("0");
    a = a < b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_test_2)
{
    BigInt a("1"), b("0");
    a = a < b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a < b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a < b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a < b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a < b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a < b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_test_1)
{
    BigInt a("0"), b("0");
    a = a > b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_test_2)
{
    BigInt a("1"), b("0");
    a = a > b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a > b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a > b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a > b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a > b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a > b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_or_equal_test_1)
{
    BigInt a("0"), b("0");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_or_equal_test_2)
{
    BigInt a("1"), b("0");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_or_equal_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_or_equal_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_or_equal_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_less_or_equal_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_less_or_equal_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a <= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_or_equal_test_1)
{
    BigInt a("0"), b("0");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_or_equal_test_2)
{
    BigInt a("1"), b("0");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_or_equal_test_3)
{
    BigInt a("-1"), b("31231234343124");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_or_equal_test_4)
{
    BigInt a("423423"), b("31231234343124");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_greater_or_equal_test_5)
{
    BigInt a("4234255555555555555553"), b("534");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_or_equal_test_6)
{
    BigInt a("67568967356264543565467347567356"), b("67568967356264543565467347567356");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_greater_or_equal_test_7)
{
    BigInt a("-67568967356264543565467347567356"), b("-67568967356264543565467347567356");
    a = a >= b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_cast_int_test_1)
{
    BigInt a("-1");
    EXPECT_EQ((int)a, -1);
}

TEST(BigIntTest, operator_cast_int_test_2)
{
    BigInt a("212121");
    EXPECT_EQ((int)a, 212121);
}

TEST(BigIntTest, operator_cast_int_test_3)
{
    BigInt a("-34654756345643785673");
    EXPECT_EQ((int)a, -643785673);
}

TEST(BigIntTest, operator_cast_string_test_1)
{
    BigInt a("-1");
    EXPECT_EQ(std :: string (a), "-1");
}

TEST(BigIntTest, operator_cast_string_test_2)
{
    BigInt a("212121");
    EXPECT_EQ(std :: string (a), "212121");
}

TEST(BigIntTest, operator_cast_string_test_3)
{
    BigInt a("-34654756345643785673");
    EXPECT_EQ(std :: string (a), "-34654756345643785673");
}

TEST(BigIntTest, operator_addition_test)
{
    BigInt a("2"), b("2");
    a = a + b;
    EXPECT_EQ(std :: string (a), "4");
}

TEST(BigIntTest, operator_subtraction_test)
{
    BigInt a("2"), b("2");
    a = a - b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_multiplication_test)
{
    BigInt a("2"), b("2");
    a = a * b;
    EXPECT_EQ(std :: string (a), "4");
}

TEST(BigIntTest, operator_division_test)
{
    BigInt a("2"), b("2");
    a = a / b;
    EXPECT_EQ(std :: string (a), "1");
}

TEST(BigIntTest, operator_xor_test)
{
    BigInt a("2"), b("2");
    a = a ^ b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_mod_test)
{
    BigInt a("2"), b("2");
    a = a % b;
    EXPECT_EQ(std :: string (a), "0");
}

TEST(BigIntTest, operator_and_test)
{
    BigInt a("2"), b("2");
    a = a & b;
    EXPECT_EQ(std :: string (a), "2");
}

TEST(BigIntTest, operator_or_test)
{
    BigInt a("2"), b("2");
    a = a | b;
    EXPECT_EQ(std :: string (a), "2");
}

TEST(BigIntTest, operator_ostream_test)
{
    BigInt a("-32423");
    std :: ostringstream strm;
    strm << a;
    EXPECT_EQ(strm.str(), "-32423");
}

TEST(BigIntTest, operator_size_test_1)
{
    BigInt a("32423");
    EXPECT_EQ(a.size(), 4);
}

TEST(BigIntTest, operator_size_test_2)
{
    BigInt a("32423543453473");
    EXPECT_EQ(a.size(), 8);
}

TEST(BigIntTest, pow_test)
{
    BigInt a("10"), b("1");
    a = pow_bigint(a, b);
    EXPECT_EQ(std :: string (a),"10");
}

TEST(BigIntTest, large_test_pow_10_10000)
{
    BigInt a("10"), b("10000");
    a = pow_bigint(a, b);
    EXPECT_EQ(std :: string (a), "10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}

TEST(BigIntTest, large_test_division_10_10000_by_2)
{
    BigInt a("10"), b("10000");
    a = pow_bigint(a, b) / BigInt(2);
    EXPECT_EQ(std :: string (a), "5000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
}