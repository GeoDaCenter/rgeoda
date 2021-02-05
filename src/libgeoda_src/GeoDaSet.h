#ifndef __GEODA_CENTER_GEODASET_H__
#define __GEODA_CENTER_GEODASET_H__

#include <string.h>

/** Old code used by LISA functions */
class GeoDaSet {
private:
    int size;
	int current;
    int* buffer;
    char* flags;
public:
	GeoDaSet(const int sz) : size(sz), current(0) {
		buffer = new int [ size ];
		flags = new char [ size ];
		memset(flags, '\x0', size);
	}
	virtual ~GeoDaSet() {
		if (buffer) delete [] buffer;
		buffer = 0;
		if (flags) delete [] flags;
		flags = 0;
		size = current = 0;
	}
    bool Belongs( const int elt) const {
		return flags[elt] != 0; }; // true if the elt belongs to the set
    void Push(const int elt) {
		// insert element in the set, if it is not yet inserted
		if (flags[elt] == 0)  {
			buffer[ current++ ] = elt;
			flags[elt] = 'i';  // check elt in
        }
    }
    int Pop() { // remove element from the set
		if (current == 0) return -1; // formerly GdaConst::EMPTY
        int rtn= buffer[ --current ];
        flags[rtn]= '\x0';   // check it out
        return rtn;
    }
    int Size() const { return current; }
	void Reset() {
		memset(flags, '\x0', size);
		current = 0;
	}
};

#endif
