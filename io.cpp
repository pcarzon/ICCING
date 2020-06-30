#include "io.h"

IO::IO()
{

}

IO::~IO()
{

}

IO::IO(const IO &original)
{

}

void IO::CopyIO(const IO &e)
{

}

IO& IO::operator= (const IO& original)
{
	CopyIO(original);
	return *this;
}
