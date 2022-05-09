#ifndef DUO_H
#define DUO_H

class duo
{
public:
    duo() {}

    duo(int _x, int _y)
    {
        xy(_x, _y);
    }

    void xy(int _x, int _y)
    {
        x(_x);
        y(_y);
    }

    void x(int _x)
    {
        if (_x < 0)
        {
            v |= (1 << 1); // sign
            v |= (1 << 0); // value
        }
        else
        {
            v &= ~(1 << 1); // sign
            if (_x == 0)
            {
                v &= ~(1 << 0); // value
            }
            else
            {
                v |= (1 << 0); // value
            }
        }
    }

    int x() const
    {
        int rs = v & (1 << 0) ? 1 : 0;
        return v & (1 << 1) ? -rs : rs;
    }

    void y(int _y)
    {
        if (_y < 0)
        {
            v |= (1 << 3); // sign
            v |= (1 << 2); // value
        }
        else
        {
            v &= ~(1 << 3); // sign
            if (_y == 0)
            {
                v &= ~(1 << 2); // value
            }
            else
            {
                v |= (1 << 2); // value
            }
        }
    }

    int y() const
    {
        int rs = v & (1 << 2) ? 1 : 0;
        return v & (1 << 3) ? -rs : rs;
    }

private:
    char v;
};

#endif // DUO_H
