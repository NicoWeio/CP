#include <cmath>
#include <iostream>
using namespace std;

float calculate_a(float x) // for large x≫1
{
    return 1 / sqrt(x) - 1 / sqrt(x + 1);
}

float calculate_b(float x) // for small x≪1
{
    return (1 - cos(x)) / sin(x);
}

float calculate_c(float x, float delta) // for small δ≪1
{
    return sin(x + delta) - sin(x);
}

int main()
{
    float result_a = calculate_a(1000000000);
    cout << "Result for a) " << result_a << endl;
    // TODO: gemeint ist |x|≪1, richtig?
    float result_b = calculate_b(0.0000000001);
    cout << "Result for b) " << result_b << endl;
    float result_c = calculate_c(1, 0.0000000001);
    cout << "Result for c) " << result_c << endl;
    return 0;
}
