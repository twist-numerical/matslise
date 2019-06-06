//
// Created by toon on 5/6/19.
//

#ifndef MATSLISE_FIND_SECTOR_H
#define MATSLISE_FIND_SECTOR_H

namespace matslise {
    template<typename Problem>
    int find_sector(const Problem *ms, typename Problem::Scalar point) {
        int a = 0, b = ms->sectorCount, c;
        while (!ms->sectors[c = a + (b - a) / 2]->contains(point)) {
            if (c == a)
                return -1;
            if (point < ms->sectors[c]->min)
                b = c;
            else
                a = c;
        }
        return c;
    }
}

#endif //MATSLISE_FIND_SECTOR_H
