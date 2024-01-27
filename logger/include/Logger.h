#pragma once

#include <memory>

#include <spdlog/spdlog.h>

namespace Logger {
    class Logger {
    public:
        static void Init();

        static std::shared_ptr<spdlog::logger>& GetLogger() noexcept;
    private:
        static std::shared_ptr<spdlog::logger> s_logger;
    };
}

#define LOG_TRACE(...) Logger::Logger::GetLogger()->trace(__VA_ARGS__)
#define LOG_INFO(...) Logger::Logger::GetLogger()->info(__VA_ARGS__)
#define LOG_DEBUG(...) Logger::Logger::GetLogger()->debug(__VA_ARGS__)
#define LOG_WARN(...) Logger::Logger::GetLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...) Logger::Logger::GetLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...) Logger::Logger::GetLogger()->critical(__VA_ARGS__)