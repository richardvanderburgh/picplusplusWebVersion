#include "Logger.h"

#include <spdlog/sinks/stdout_color_sinks.h>

namespace Logger {
	std::shared_ptr<spdlog::logger> Logger::s_logger;

	void Logger::Init()
	{
		spdlog::set_pattern("*** [%n] | [%l] | [%H:%M:%S] [thread %t] | %v ***");
		s_logger = spdlog::stdout_color_mt("Core Logger");
		s_logger->set_level(spdlog::level::trace);
	}

	std::shared_ptr<spdlog::logger>& Logger::GetLogger() noexcept
	{
		return s_logger;
	}
}

