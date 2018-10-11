/*
 * LtrTe.cpp
 *
 *  Created on: Dec 26, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "LtrTe.h"
#include "BackwardTr.h"

#include "../exception/InvalidOperationException.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"
#include "../utility/ITail.h"
#include "../utility/Tail.h"
#include "../utility/EmptyTail.h"
#include "../utility/EmptyTSD.h"

using namespace utility;
using namespace exception;

namespace tr {

LtrTe::LtrTe(BackwardTr * ltrIn, ITSD * tsdIn, ITail * pptIn) {
	initializer(ltrIn, tsdIn, pptIn);
}

LtrTe::LtrTe(LtrTe& copy) {
	initializer(copy.getLtr(), copy.getTsd(), copy.getPpt());
}

void LtrTe::initializer(BackwardTr * ltrIn, ITSD * tsdIn, ITail * pptIn) {
	ltr = new BackwardTr(*ltrIn);

	setTsd(tsdIn);
	setPpt(pptIn);

	if (s > e) {
		string msg("The start of the TE must be <= its end. ");
		msg.append("The start is: ");
		msg.append(Util::int2string(s));
		msg.append(" The end is: ");
		msg.append(Util::int2string(e));
		msg.append(".");
		throw InvalidStateException(msg);
	}

	if (tsd != EmptyTSD::getInstance()
			&& !(ltr->getStart() > s && ltr->getEnd() < e)) {
		string msg("The LTR must be within the TSD. ");
		msg.append("The start and the end of the TSD are: ");
		msg.append(Util::int2string(s));
		msg.append(":");
		msg.append(Util::int2string(e));
		msg.append(". The start and the end of the LTR are: ");
		msg.append(Util::int2string(ltr->getStart()));
		msg.append(":");
		msg.append(Util::int2string(ltr->getEnd()));
		msg.append(".");
		throw InvalidStateException(msg);
	}
}

LtrTe::~LtrTe() {
	delete ltr;
	if (tsd != EmptyTSD::getInstance()) {
		delete tsd;
	}

	if (ppt != EmptyTail::getInstance()) {
		delete ppt;
	}
}

int LtrTe::getStart() const {
	return s;
}

int LtrTe::getEnd() const {
	return e;
}

int LtrTe::getLength() {
	return e - s + 1;
}

BackwardTr * LtrTe::getLtr() {
	return ltr;
}

ITail * LtrTe::getPpt() {
	return ppt;
}

ITSD * LtrTe::getTsd() {
	return tsd;
}
//11/9/17 changed to print in .bed format rather than .coor
string LtrTe::toString(string header) {
	string msg("");
	msg.append(header);
	msg.append("	");
	msg.append(toStringHelper());
	return msg;
}

string LtrTe::toString() {
	string msg("LTR_TE ");
	msg.append(toStringHelper());
	return msg;
}

string LtrTe::toStringHelper() {
	string msg(Util::int2string(s));
	msg.append("	");
	msg.append(Util::int2string(e));
	/*msg.append(" ");
	msg.append(ltr->toString());
	msg.append(" ");
	msg.append(tsd->toString());
	msg.append(" ");
	msg.append(ppt->toString());*/
	return msg;
}

void LtrTe::setStart(int start) {
	string msg("Setting the start of an instance of LTR TE is not allowed.");
	throw InvalidOperationException(msg);
}

void LtrTe::setEnd(int end) {
	string msg("Setting the end of an instance of LTR TE is not allowed.");
	throw InvalidOperationException(msg);
}

void LtrTe::setPpt(ITail* pptIn) {
	if (pptIn == EmptyTail::getInstance()) {
		ppt = pptIn;
	} else {
		ppt = new Tail(*pptIn);
	}
}

void LtrTe::setTsd(ITSD* tsdIn) {
	if (tsdIn == EmptyTSD::getInstance()) {
		tsd = tsdIn;
		s = ltr->getStart();
		e = ltr->getEnd();
		
	} else {
		tsd = new TSD(*tsdIn);
		s = tsd->getLtTsd()->getStart();
		e = tsd->getRtTsd()->getEnd();
	}
}

bool LtrTe::lessThan(LtrTe* a, LtrTe* b) {
	return a->getStart() < b->getStart();
}

} /* namespace tr */
