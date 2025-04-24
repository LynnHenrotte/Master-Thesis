from typing import Final
from dataclasses import dataclass

@dataclass(frozen=True)
class APIConfig:
    """API Configuration settings"""
    DEFAULT_HOST: Final[str] = "www.pharmvar.org/api-service"
    DEFAULT_VERSION: Final[str] = "0.1"
    DEFAULT_SSL_VERIFY: Final[bool] = False
    DEFAULT_LOGGER = None

@dataclass(frozen=True)
class LogConfig:
    """Logging Configuration settings"""
    DEFAULT_LOG_LEVEL: Final[str] = "INFO"
    DEFAULT_LOG_FORMAT: Final[str] = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Package metadata
__version__: Final[str] = "0.1.0"
__author__: Final[str] = "Your Name"
__author_email__: Final[str] = "your.email@example.com"
__description__: Final[str] = "A Python client for the PharmVar API"

# HTTP settings
TIMEOUT: Final[int] = 30  # seconds
MAX_RETRIES: Final[int] = 3
RETRY_BACKOFF_FACTOR: Final[float] = 0.5

# Rate limiting
RATE_LIMIT: Final[int] = 100  # requests per minute

# Cache settings
CACHE_TTL: Final[int] = 3600  # seconds
CACHE_ENABLED: Final[bool] = True