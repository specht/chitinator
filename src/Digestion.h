/*
Copyright (c) 2007-2008 Michael Specht

This file is part of Chitinator.

Chitinator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Chitinator.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <QtCore>
#include "RefPtr.h"
#include "Polymer.h"
#include "Enzyme.h"


class k_Digestion
{
public:
	k_Digestion();
	virtual ~k_Digestion();
	
	virtual QList<RefPtr<k_Polymer> > digest(QList<RefPtr<k_Polymer> > ak_Polymers, k_Enzyme ak_Enzyme);
};
