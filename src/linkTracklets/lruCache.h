// jmyers: took this from http://datahack.se/?p=103
// when done with testing replace it with somebody's real implementation!

#ifndef _LRUCACHE_HPP_
#define _LRUCACHE_HPP_

#include <list>
#include <map>

template<typename Key, typename Value>
class LRUCache
{
public:
     /* Shortcuts */

     typedef std::list< std::pair< Key, Value> > List;
     typedef typename List::iterator ListIter;
     typedef typename List::const_iterator ListIterConst;
     typedef std::map< Key, ListIter > Index;
     typedef typename Index::iterator IndexIter;
     
     /* Constructor */
     
     LRUCache(size_t size)  {
	  m_size = size;
     }

     ~LRUCache() {
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

     bool empty() const {
	  return m_list.empty();
     }

     size_t size() const {
	  return m_list.size();
     }

     void set_max_size(size_t size) {
	  m_size = size;
     }

     size_t max_size() const {
	  return m_size;
     }

     /* Modifiers */

     void insert(const Key& k, const Value& v) {
	  IndexIter _i = m_index.find(k);

	  // update current item
	  if (_i != m_index.end()) {
	       (*_i->second).second = v;
	       _touch(_i);
	       return;
	  }

	  // insert new item
	  m_index.insert(
	       std::make_pair(
		    k,
		    m_list.insert(
			 m_list.begin(),
			 std::make_pair(k, v)
			 )
		    )
	       );

	  // truncate if the list is too big
	  if (m_size && m_list.size() > m_size) {
	       ListIter _l = m_list.end();
	       --_l;
	       remove(_l->first);
	  }
     }

     void remove(const Key& k) {
	  IndexIter _i = m_index.find(k);
	  if (_i != m_index.end()) {
	       m_list.erase(_i->second);
	       m_index.erase(_i);
	  }
     }

     void clear() {
	  m_index.clear();
	  m_list.clear();
     }

     /* Operations */

     bool find(const Key& k, Value& v, bool touch = true) {
	  IndexIter _i = m_index.find(k);
	  if (_i != m_index.end()) {
	       v = (*_i->second).second;
	       if (touch) _touch(_i);
	       return true;
	  }
	  return false;
     }

protected:
     size_t m_size;
     List m_list;
     Index m_index;

     // move item to head
     void _touch(IndexIter& _i) {
	  m_list.splice(m_list.begin(), m_list, _i->second);
	  _i->second = m_list.begin();
     }
};

#endif
