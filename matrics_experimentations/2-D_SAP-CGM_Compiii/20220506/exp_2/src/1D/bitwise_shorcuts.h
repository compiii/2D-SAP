#ifndef BITWISE_SHORCUTS_H
#define BITWISE_SHORCUTS_H

// USE WITH CAUTION: vlong offsets must be consistent

#include <functional>

struct vlong
{
    long sign = 0;
    long value = 0;

    bool operator==(const vlong &o) const
    {
        return (sign == o.sign) && (value == o.value);
    }
};

namespace std
{
    template <>
    class hash<vlong>
    {
    public:
        size_t operator()(const vlong &offsets) const
        {
            return std::hash<long>{}(offsets.sign + offsets.value);
        }
    };
} // namespace std

inline int bit_get(vlong offsets, int pos)
{
    if (offsets.sign & (1 << pos))
    {
        return -1;
    }

    return (offsets.value & (1 << pos)) ? 1 : 0;
}

inline void bit_set(vlong &offsets, int pos, int val)
{
    if (val < 0)
    {
        offsets.sign |= (1 << pos);
        offsets.value |= (1 << pos);
        return;
    }

    // offsets.sign  &= ~(1 << pos);
    // offsets.value &= ~(1 << pos);

    if (val > 0)
    {
        offsets.value |= (1 << pos);
    }
}

#endif // BITWISE_SHORCUTS_H
