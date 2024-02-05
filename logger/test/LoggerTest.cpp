#include <gtest/gtest.h>

#include <Logger.h>

TEST(LoggerTest, create_logger) {
	EXPECT_TRUE(Logger::Logger::GetLogger() == nullptr);
	Logger::Logger::Init();
	EXPECT_FALSE(Logger::Logger::GetLogger() == nullptr);
}