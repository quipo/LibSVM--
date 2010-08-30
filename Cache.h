/* 
 * File:   Cache.h
 *
 * Created on July 3, 2010, 3:30 PM
 */

#ifndef _CACHE_HPP_
#define _CACHE_HPP_

#include <list>
#include <map>

template<typename Key, typename Value>
class Cache {
public:
    /* Shortcuts */

    typedef std::list< std::pair< Key, Value> > List;
    typedef typename List::iterator ListIter;
    typedef typename List::const_iterator ListIterConst;
    typedef std::map< Key, ListIter > Index;
    typedef typename Index::iterator IndexIter;

    /* Constructor */

    Cache() : m_size(0) {
    }

    Cache(size_t size) : m_size(size) {
    }

    virtual ~Cache() {
    }

    /* Iterators */

    ListIterConst begin() const {
        return m_list.begin();
    }

    ListIterConst rbegin() const {
        return m_list.rbegin();
    }

    ListIterConst end() const {
        return m_list.end();
    }

    ListIterConst rend() const {
        return m_list.rend();
    }

    /* Capacity */

    virtual bool empty() const = 0;

    virtual size_t size() const = 0;

    virtual void set_max_size(size_t size) = 0;

    virtual size_t max_size() const = 0;

    /* Modifiers */

    virtual void insert(const Key& k, const Value& v) = 0;

    virtual void remove(const Key& k) = 0;

    virtual void clear() = 0;

    /* Operations */

    virtual bool find(const Key& k, Value& v, bool touch = true) = 0;

    /* Lorenzo */

    virtual void swap_index(const Key i, const Key j) = 0;

protected:
    size_t m_size;
    List m_list;
    Index m_index;
};

#endif /* _CACHE_HPP_ */
