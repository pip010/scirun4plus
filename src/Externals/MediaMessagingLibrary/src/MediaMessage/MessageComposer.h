#ifndef MESSAGE_COMPOSER_H
#define MESSAGE_COMPOSER_H

#include <boost/noncopyable.hpp>
#include <boost/thread/tss.hpp>

#include "Message.h"

namespace BioMesh3d
{
	class MessageComposer : private boost::noncopyable
	{
	public:
		MessageComposer() {}
		~MessageComposer() {}
	protected:
	};
}

#endif