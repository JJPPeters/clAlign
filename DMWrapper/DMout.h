#pragma once

#include "stdafx.h"

static struct DMOutResult
{
	DMOutResult& operator<<(std::string value)
	{
		DM::Result(value);
		return *this;
	}

	template<typename T>
	DMOutResult& operator<<(T value)
	{
		std::stringstream ss;
		ss << value;
		DM::Result(ss.str());
		return *this;
	}
} DMresult;

// Very, very bad way of defining an endline
static std::string DMendl = "\n";