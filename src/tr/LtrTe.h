/*
 * LtrTe.h
 *
 *  Created on: Dec 26, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef LTRTE_H_
#define LTRTE_H_

#include "BackwardTr.h"
#include "../utility/ITSD.h"
#include "../utility/TSD.h"
#include "../utility/ITail.h"

using namespace utility;

namespace tr {

class LtrTe: public ILocation {
private:
	BackwardTr * ltr;
	ITSD * tsd;
	ITail * ppt;
	int s;
	int e;
	void initializer(BackwardTr *, ITSD *, ITail *);
	string toStringHelper();

public:
	LtrTe(BackwardTr *, ITSD *, ITail *);
	LtrTe(LtrTe &);
	virtual ~LtrTe();

	virtual int getEnd() const;
	virtual int getStart() const;
	virtual void setEnd(int);
	virtual void setStart(int);
	virtual int getLength();
	virtual string toString();
	string toString(string);

	BackwardTr * getLtr();
	ITail * getPpt();
	ITSD * getTsd();

	void setPpt(ITail *);
	void setTsd(ITSD *);

	static bool lessThan(LtrTe* a, LtrTe* b);

};

} /* namespace tr */
#endif /* LTRTE_H_ */
