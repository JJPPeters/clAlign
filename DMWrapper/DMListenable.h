#pragma once

// used so listeners can take any object and call it's process method.
class DMListenable abstract
{
public:
	virtual void Process()=0;
};