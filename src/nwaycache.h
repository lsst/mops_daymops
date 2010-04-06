// Copyright (c) 2004 Ivan Skytte Jørgensen
//
// This software is provided 'as-is', without any express or implied warranty.
// In no event will the authors be held liable for any damages arising from the
// use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution.


#ifndef _NWAYCACHE_H_
#define _NWAYCACHE_H_

// An N-Way set associative cache
// Features:
//   Linear memory
//   Automatic purging of old items
//   Supports both caching presence and non-presence of items

#include <stddef.h>
#include <time.h>
#include <new>

namespace NWayCache {

//Forward declarations ad nauseam
template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
class Cache;
template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
class SlotReference;
template <class K, class T>
struct Slot;
template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
class CacheIterator;

//Default no-op cache and bucket locks
class NullLockObject {
	NullLockObject(const NullLockObject&);
	NullLockObject& operator=(const NullLockObject&);
public:
	NullLockObject() {}
//Not needed becausee the Locker is specialized for NullLockObject:
//	void lock() {}
//	void unlock() {}
};


//A simple scoped lock
template<class LockObject>
class Locker {
	LockObject *lo;
	Locker(const Locker&); //don't copy
	Locker& operator=(const Locker&); //don't copy
public:
	Locker() : lo(0) {}
	Locker(LockObject &lo_) : lo(&lo_) { lo->lock(); }
	~Locker() { release(); }
	void lock(LockObject *lo_) { lo=lo_; lo->lock(); }
	void release() { if(lo) { lo->unlock(); lo=0; } }
	void swap(Locker<LockObject> &l) {
		LockObject *tmp = lo;
		lo = l.lo;
		l.lo = tmp;
	}
};

//Specialized no-op locker.
template<>
class Locker<NullLockObject> {
public:
	Locker()  {}
	Locker(NullLockObject &) { }
	~Locker() {  }
	void lock(NullLockObject *) { }
	void release() { }
	void swap(Locker<NullLockObject> &) { }
};



template <class K, class T>
struct Slot {
	bool	in_use;
	bool	available;
	time_t	expire_timestamp;
	K	key;
	T	data;
	
	void clear() { if(in_use) key.~K(); if(available) data.~T(); in_use=available=false; }
};

template <class K, class T, size_t N, class LockObject>
struct Bucket {
	typedef Slot<K,T> Slot_t;
	LockObject lock_object;
	union {
		char	slot_buf[sizeof(Slot_t)*N];
		//needed to force alignment of slot_buf
		void	*ptr_align_;
		double	double_align_;
		unsigned long ul_align_;
	};
	size_t	lru[N];
	Slot_t	&slot(size_t slot_no) { return ((Slot_t*)slot_buf)[slot_no]; }
	typedef Slot_t *iterator;
	iterator begin() { return (iterator)slot_buf; }
	iterator end() { return begin()+N; }
	
	Bucket() {
		for(size_t i=0; i<N; i++) {
			slot(i).in_use = false;
			lru[i] = i;
		}
	}
	~Bucket() {
		for(size_t i=0; i<N; i++)
			slot(i).clear();
	}
	Slot<K,T> *find(const K &key, time_t now);
	void update_lru(Slot_t *s);
	Slot<K,T> *selectVictimSlot() { return &slot(lru[N-1]); }
	size_t size() const { size_t sz; for(size_t i=0; i<N; i++) if(slot(i).in_use) sz++; return sz; }
	void clear() { for(size_t i=0; i<N; i++) slot(i).clear(); }
};


struct NotAvailable_t {
};
static const NotAvailable_t NotAvailable = {};


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
class SlotReference {
	typedef Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject> Cache_t;
	Cache_t *cache;
	typedef Bucket<K,T,N,BucketLockObject> Bucket_t;
	Bucket_t *bucket;
	typedef Slot<K,T> Slot_t;
	Slot_t *slot;
	
	Locker<CacheLockObject> cache_locker;
	Locker<BucketLockObject> bucket_locker;
	friend class Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>;
	void unset() {
		cache=0;
		bucket=0;
		slot=0;
		cache_locker.release();
		bucket_locker.release();
	}
	void set(Cache_t *cache_,
	         Bucket_t *bucket_,
	         Slot_t *slot_,
		 Locker<CacheLockObject> &cl,
		 Locker<BucketLockObject> &bl)
	{
		unset();
		cache = cache_;
		bucket = bucket_;
		slot = slot_;
		cache_locker.swap(cl);
		bucket_locker.swap(bl);
	}
public:
	SlotReference() : cache(0), bucket(0), slot(0) {}
	operator bool() const { return slot!=0; }
	bool available() const { return slot->available; }
	T *operator->() { return &slot->data; }
	T& operator*() { return slot->data; }
	const K& key() const { return slot->key; }
	const T& data() const { return slot->data; }
	bool operator==(const NotAvailable_t &) { return !slot->available; }
	bool operator!=(const NotAvailable_t &) { return slot->available; }
	void release() { unset(); }
	void update(const T& data);
	void update(const NotAvailable_t&);
	void update();
	void erase();
};

template<class K>
struct GenericHasher {
static unsigned long hash(const K& k) { return (unsigned long)k; }
};


template <class K, class T, size_t N, class Hasher=GenericHasher<K>, class CacheLockObject=NullLockObject, class BucketLockObject=NullLockObject>
class Cache {
public:
	typedef SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject> SlotReference_t;
	typedef CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject> iterator;
private:
	typedef Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject> Self_t;
	typedef Bucket<K,T,N,BucketLockObject> Bucket_t;
	typedef typename Bucket_t::Slot_t Slot_t;
	typedef Locker<CacheLockObject> CacheLocker_t;
	typedef Locker<BucketLockObject> BucketLocker_t;
	friend class Slot<K,T>;
	friend class SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject>;
	
	CacheLockObject lock_object;
	const size_t num_buckets;
	Bucket_t *buckets;
	const time_t max_item_age;
	
	//Statistics are not guaranteed to be consistent unless cachelockobject synchronizes
	unsigned long num_positive_hits;
	unsigned long num_negative_hits;
	unsigned long num_misses;
	
	Cache(const Self_t&); //don't copy
	Self_t& operator=(const Self_t&); //don't copy
	
	bool insert_(const K &key, const T *data);
	bool update_(const K &key, const T *data);
	void update_(SlotReference_t &sr, const T *data);
	void crupdate_(const K &key, const T *data);
	
public:
	Cache(size_t num_buckets_, time_t max_item_age_)
	  : lock_object(),
	    num_buckets(num_buckets_),
	    buckets(new Bucket_t[num_buckets]),
	    max_item_age(max_item_age_),
	    num_positive_hits(0),
	    num_negative_hits(0),
	    num_misses(0)
	  { }
	~Cache() { delete[] buckets; }
	
	bool find(SlotReference_t &sr, const K &key);
	bool exists(const K &key);
	void erase(const K &key);
	void erase(SlotReference_t &sr);
	bool insert(const K &key, const T &data) { return insert_(key,&data); }
	bool insert(const K &key, const NotAvailable_t&) { return insert_(key,0); }
	bool update(const K &key, const T &data) { return update_(key,&data); }
	bool update(const K &key, const NotAvailable_t&) { return update_(key,0); }
	void update(SlotReference_t &sr, const T &data) { update_(sr,&data); }
	void update(SlotReference_t &sr, const NotAvailable_t&) { update_(sr,0); }
	void crupdate(const K &key, const NotAvailable_t&) { crupdate_(key,0); }
	void crupdate(const K &key, const T &data) { crupdate_(key,&data); }
	
	void queryStatistics(unsigned long &num_positive_hits_,
	                     unsigned long &num_negative_hits_,
			     unsigned long &num_misses_,
			     bool reset=false)
	{
		num_positive_hits_ = num_positive_hits;
		num_negative_hits_ = num_negative_hits;
		num_misses_ = num_misses;
		if(reset) {
			num_positive_hits=0;
			num_negative_hits=0;
			num_misses=0;
		}
	}
	
	size_t queryNumBuckets() const { return num_buckets; }
	
	//STL helpers
	size_t max_size() const { return num_buckets*N; }
	size_t capacity() const { return num_buckets*N; }
	bool empty() const { return size()==0; }
	size_t size() const;
	void clear();
	
	iterator begin();
	iterator end();
	iterator erase(iterator &ci);
	void advance_(iterator *i);
};


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
class CacheIterator {
public:
	typedef CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject> Self_t;
	struct pair {
		const K &first;
		T &second;
		pair(const K& key, T &data) : first(key), second(data) {}
	};
private:
	char pair_buf[sizeof(pair)];
	pair& the_pair;
	typedef Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject> Cache_t;
	Cache_t *cache;
	typedef Bucket<K,T,N,BucketLockObject> Bucket_t;
	Bucket_t *bucket;
	typedef Slot<K,T> Slot_t;
	typename Bucket_t::iterator slot;
	
	void set(Cache_t *cache_, Bucket_t *bucket_, Slot_t *slot_);
	friend class Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>;
	CacheIterator(Cache_t *cache_, Bucket_t *bucket_, Slot_t *slot_)
	  : the_pair(*(pair*)pair_buf)
	  { set(cache_,bucket_,slot_); }
public:
	CacheIterator()
	  : the_pair(*(pair*)pair_buf)
	  {}
	CacheIterator(const Self_t& ci)
	  : the_pair(*(pair*)pair_buf)
	  { set(ci.cache,ci.bucket,ci.slot); }
	Self_t& operator=(const Self_t& ci)
	  { set(ci.cache,ci.bucket,ci.slot); return *this; }
	bool operator==(const Self_t &ci) const { return slot==ci.slot; }
	bool operator!=(const Self_t &ci) const { return slot!=ci.slot; }
	
	bool available() const { return slot->available; }
	pair* operator->() { return &the_pair; }
	pair& operator*() { return the_pair; }
	Self_t& operator++() { cache->advance_(this); return *this; }
	Self_t operator++(int) { Self_t tmp=*this; ++*this; return tmp; }
};


//-----------------------------------------------------------------------------


template <class K, class T, size_t N, class BucketLockObject>
inline Slot<K,T> *Bucket<K,T,N,BucketLockObject>::find(const K &key, time_t now) {
	for(size_t i=0; i<N; i++) {
		Slot_t *s = &slot(i);
		if(s->in_use && s->key == key) {
			if(s->expire_timestamp<now) {
				s->clear();
				return false;
			}
			return s;
		}
	}
	return 0;
}

inline void update_lru_(size_t lru[], size_t i) {
	if(i!=lru[0]) {
		size_t j=0;
		while(lru[j]!=i)
			j++;
		while(j!=0) {
			lru[j]=lru[j-1];
			j--;
		}
		lru[0]=i;
	}
}

template <class K, class T, size_t N, class BucketLockObject>
inline void Bucket<K,T,N,BucketLockObject>::update_lru(Slot_t *s) {
	size_t i = (size_t)(s-&slot(0));
	update_lru_(lru,i);
}



template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject>::update() {
	if(slot) slot->expire_timestamp=time(0)+cache->max_item_age;
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject>::update(const T& data) {
	cache->update(*this,data);
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject>::update(const NotAvailable_t&) {
	cache->update(*this,NotAvailable);
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void SlotReference<K,T,N,Hasher,CacheLockObject,BucketLockObject>::erase() {
	if(slot) cache->erase(*this);
}



template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline bool Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::find(SlotReference_t &sr, const K &key) {
	time_t now = time(0);
	unsigned long hash_value = Hasher::hash(key);
	
	size_t bucket_idx = hash_value%num_buckets;
	Bucket_t *bucket = &buckets[bucket_idx];
	CacheLocker_t cl(lock_object);
	BucketLocker_t bl(bucket->lock_object);
	Slot_t *slot = bucket->find(key,now);
	if(!slot) {
		num_misses++;
		return false;
	}
	if(slot->available)
		num_positive_hits++;
	else
		num_negative_hits++;
	
	sr.set(this,bucket,slot,cl,bl);
	return true;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline bool Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::exists(const K &key) {
	time_t now = time(0);
	unsigned long hash_value = Hasher::hash(key);
	
	size_t bucket_idx = hash_value%num_buckets;
	Bucket_t *bucket = &buckets[bucket_idx];
	CacheLocker_t cl(lock_object);
	BucketLocker_t bl(bucket->lock_object);
	Slot_t *slot = bucket->find(key,now);
	
	if(!slot) {
		num_misses++;
		return false;
	}
	
	if(slot->available)
		num_positive_hits++;
	else
		num_negative_hits++;
	return true;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::erase(const K &key) {
	SlotReference_t sr;
	if(find(sr,key)) erase(sr);
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::erase(SlotReference_t &sr) {
	if(!sr) return;
	sr.slot->clear();
	sr.release();
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline bool Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::insert_(const K &key, const T *data) {
	time_t now = time(0);
	unsigned long hash_value = Hasher::hash(key);
	
	size_t bucket_idx = hash_value%num_buckets;
	Bucket_t *bucket = &buckets[bucket_idx];
	CacheLocker_t cl(lock_object);
	BucketLocker_t bl(bucket->lock_object);
	Slot_t *slot = bucket->find(key,now);
	if(slot)
		return false; //will not overwrite
	slot = bucket->selectVictimSlot();
	if(slot->in_use)
		slot->key = key;
	else {
		new(&slot->key) K(key);
		slot->in_use = true;
	}
	if(data) {
		if(slot->available)
			slot->data = *data;
		else {
			new(&slot->data) T(*data);
			slot->available=true;
		}
	} else {
		if(slot->available) {
			slot->data.~T();
			slot->available=false;
		}
	}
	slot->expire_timestamp = now + max_item_age;
	bucket->update_lru(slot);
	return true;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline bool Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::update_(const K &key, const T *data) {
	SlotReference_t sr;
	if(!find(sr,key)) return false;
	update_(sr,data);
	return true;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::update_(SlotReference_t &sr, const T *data) {
	if(!sr) return;
	time_t now = time(0);
	
	if(data) {
		if(sr.slot->available)
			sr.slot->data = *data;
		else {
			new(&sr.slot->data) T(*data);
			sr.slot->available=true;
		}
	} else {
		if(sr.slot->available) {
			sr.slot->data.~T();
			sr.slot->available=false;
		}
	}
	sr.slot->expire_timestamp = now + max_item_age;
	sr.bucket->update_lru(sr.slot);
	return;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::crupdate_(const K &key, const T *data) {
	time_t now = time(0);
	unsigned long hash_value = Hasher::hash(key);
	
	size_t bucket_idx = hash_value%num_buckets;
	Bucket_t *bucket = &buckets[bucket_idx];
	CacheLocker_t cl(lock_object);
	BucketLocker_t bl(bucket->lock_object);
	Slot_t *slot = bucket->find(key,now);
	if(!slot) {
		//insert
		slot = bucket->selectVictimSlot();
		if(slot->in_use)
			slot->key = key;
		else {
			new(&slot->key) K(key);
			slot->in_use = true;
		}
	}
	if(data) {
		if(slot->available)
			slot->data = *data;
		else {
			new(&slot->data) T(*data);
			slot->available=true;
		}
	} else {
		if(slot->available) {
			slot->data.~T();
			slot->available=false;
		}
	}
	slot->expire_timestamp = now + max_item_age;
	bucket->update_lru(slot);
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline size_t Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::size() const {
	size_t sz=0;
	CacheLocker_t cl(lock_object);
	for(size_t bucket_idx=0; bucket_idx<size; bucket_idx++) {
		Bucket_t *bucket = &buckets[bucket_idx];
		BucketLocker_t bl(bucket->lock_object);
		sz += bucket->size();
	}
	return sz;
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::clear() {
	CacheLocker_t cl(lock_object);
	for(size_t bucket_idx=0; bucket_idx<size; bucket_idx++) {
		Bucket_t *bucket = &buckets[bucket_idx];
		BucketLocker_t bl(bucket->lock_object);
		bucket->clear();
	}
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject> Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::begin() {
	for(Bucket_t *b=buckets; b<buckets+num_buckets; b++) {
		for(size_t i=0; i<N; i++) {
			Slot_t *slot = &b->slot(i);
			if(slot->in_use)
				return iterator(this,b,slot);
		}
	}
	return end();
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject> Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::end() {
	iterator i;
	i.cache=this;
	i.bucket=buckets+num_buckets-1;
	i.slot=i.bucket->end();
	return i;
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject> Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::erase(iterator &ci) {
	iterator tmp=ci;
	erase(ci->key);
	advance_(&ci);
	return ++tmp;
}

template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void Cache<K,T,N,Hasher,CacheLockObject,BucketLockObject>::advance_(iterator *i) {
	Bucket_t *bucket = i->bucket;
	Slot_t *slot=i->slot;
	slot++;
	while(slot!=bucket->end()) {
		if(slot->in_use) {
			i->set(this,bucket,slot);
			return;
		}
		slot++;
	}
	bucket++;
	for(;bucket<buckets+num_buckets; bucket++) {
		slot=bucket->begin();
		while(slot!=bucket->end()) {
			if(slot->in_use) {
				i->set(this,bucket,slot);
				return;
			}
			slot++;
		}
	}
	*i = end();
}


template <class K, class T, size_t N, class Hasher, class CacheLockObject, class BucketLockObject>
inline void CacheIterator<K,T,N,Hasher,CacheLockObject,BucketLockObject>::set(Cache_t *cache_, Bucket_t *bucket_, Slot_t *slot_) {
	cache=cache_;
	bucket=bucket_;
	slot=slot_;
	if(slot && slot!=bucket->end())
		new(pair_buf) pair(slot->key,slot->data);
}


} //namespace

#endif


#ifdef PTHREAD_MUTEX_INITIALIZER
#ifndef _NWAYCACHE_PTHREAD_MUTEX_LOCK_OBJECT_
#define _NWAYCACHE_PTHREAD_MUTEX_LOCK_OBJECT_
namespace NWayCache {
class MutexLockObject {
	pthread_mutex_t mtx;
	MutexLockObject(const MutexLockObject&); //don-t copy
	MutexLockObject& operator=(const MutexLockObject&); //dont copy.
public:
	MutexLockObject() { pthread_mutex_init(&mtx,0); }
	~MutexLockObject() { pthread_mutex_destroy(&mtx); }
	void lock() { pthread_mutex_lock(&mtx); }
	void unlock() { pthread_mutex_unlock(&mtx); }
};
} //namespace
#endif
#endif


#ifdef NWAYCACHE_SUPPORT_STRING
#ifndef _NWAYCACHE_STRING_HASHER_
#define _NWAYCACHE_STRING_HASHER_
#include <string>
namespace NWayCache {
template<>
struct GenericHasher<std::string> {
	static unsigned long hash(const std::string &s) {
		unsigned long v=0;
		for(std::string::const_iterator i=s.begin();
		    i!=s.end();
		    ++i)
		{
			char c = *i;
			v = v*5 + (unsigned long)c;
		}
		return v;
	}
};
} //namespace
#endif
#endif
