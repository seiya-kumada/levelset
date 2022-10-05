#ifndef LevelSetMethod2_Utilities_h
#define LevelSetMethod2_Utilities_h

namespace lsm
{
        constexpr int power(int x, int y)
        {
                return y == 1 ? x : x * power(x, y - 1);
        }
} // namespace lsm
#endif //  LevelSetMethod2_Utilities_h
