/* 
 * File:   LRUCache.cpp
 * Author: lorenzo
 * 
 * Created on July 3, 2010, 3:30 PM
 */

#include "LRUCache.h"

/*
It should be implemented like this..

#include "lru.hpp"

#include <string>
#include <stdio.h>

int main() {
    LRUCache<std::string, int> lru;
    lru.set_max_size(2);
    lru.insert("foo", 1);
    lru.insert("bar", 2);
    lru.insert("foo", 3);
    lru.insert("moo", 4);

    int f;
    if (lru.find("foo", f))
        printf("found foo %d\n", f);
    if (lru.find("bar", f, false))
        printf("found bar %d\n", f);
    if (lru.find("moo", f, false))
        printf("found moo %d\n", f);

    for (LRUCache<std::string, int>::ListIterConst l = lru.begin();
        l != lru.end(); l++) {
        printf("item %d\n", l->second);
    }

    return 0;
}

*/
