import logging
import requests
import requests.packages
from typing import List, Dict
from exceptions.exceptions import PharmVarApiException, NoDataFoundError
from models import Result
from json import JSONDecodeError
from config import APIConfig


class RestAdapter:
    def __init__(self, hostname: str = APIConfig.DEFAULT_HOST , api_key: str = "", version: str = APIConfig.DEFAULT_VERSION, logger: logging.Logger = APIConfig.DEFAULT_LOGGER):
        """
        :param hostname: The hostname of the API server: e.g. www.pharmvar.org/api-service
        :param api_key (optional): The API key to use for authentication
        :param version: The version of the API to use: currently only "0.1" is supported
        :param logger (optional): pass your logger here to use it, otherwise a new logger will be created
        """

        self.url = f"https://{hostname}/"
        self._api_key = api_key
        self._logger = logger or logging.getLogger(__name__)
        self._version = version

    def _do(self, http_method: str, endpoint: str, params: Dict = None, data: Dict = None, headers: Dict = {"Accept": "*/*"}) -> Result:
        """
        Execute HTTP request with logging and error handling
        
        Args:
            http_method (str): HTTP method (GET, POST, etc.)
            endpoint (str): API endpoint
            params (Dict, optional): Query parameters
            data (Dict, optional): Request body data
            headers (Dict, optional): Request headers
                
        Returns:
            Result: Response data wrapped in Result object
                
        Raises:
            PharmVarApiException: If request fails or response is invalid
            NoDataFoundError: If no data is returned for valid request
        """
        full_url = f"{self.url}{endpoint}"
        log_line_pre = f"method={http_method}, url={full_url}, params={params}"
        
        try:
            self._logger.debug(msg=log_line_pre)   
            response = requests.request(
                method=http_method,
                url=full_url,
                headers=headers,
                params=params,
                data=data
            )
        except requests.exceptions.RequestException as e:
            self._logger.error(msg=str(e))
            raise PharmVarApiException("Request failed") from e

        # Check content type of response
        content_type = response.headers.get('Content-Type', '')
        
        # Try to parse response data based on content type
        try:
            if 'application/json' in content_type or headers.get('Accept') == 'application/json':
                data_out = response.json()
            elif 'text/plain' in content_type or headers.get('Accept') == 'text/plain':
                data_out = response.text
                # Handle case where error response is JSON even with text/plain
                if response.status_code >= 400:
                    try:
                        data_out = response.json()
                    except:
                        pass
            else:
                # For */* or unknown content types, try JSON first then fall back to text
                try:
                    data_out = response.json()
                except (JSONDecodeError, ValueError):
                    data_out = response.text
        except (JSONDecodeError, ValueError) as e:
            log_msg = f"{log_line_pre}, success=False, status_code=None, message={str(e)}"
            self._logger.error(msg=log_msg)
            raise PharmVarApiException("Failed to parse response data") from e

        # Check for success and handle errors
        is_success = 200 <= response.status_code <= 299
        log_msg = f"{log_line_pre}, success={is_success}, status_code={response.status_code}, message={response.reason}"
        
        if is_success:
            self._logger.debug(msg=log_msg)
            result = Result(status_code=response.status_code, message=response.reason, data=data_out)
            if not result.data and isinstance(data_out, (list, dict)):
                raise NoDataFoundError(f"No data found for endpoint: {endpoint}")
            return result

        # Handle error response
        self._logger.error(msg=log_msg)
        
        # Extract error message from response
        error_message = None
        if isinstance(data_out, dict):
            error_message = data_out.get("errorMessage")
        if not error_message:
            error_message = response.reason
            
        # Special handling for 404 errors
        if response.status_code == 404 or (isinstance(data_out, dict) and data_out.get("errorCode") == 404):
            raise NoDataFoundError(error_message)
        
    def get(self, endpoint: str, params: Dict = None) -> Result:
        """
        Perform an HTTP GET request to the API
        :param endpoint: The API endpoint to call
        :param params: A dictionary of parameters to pass to the API
        :return: A Result object containing the status code, message, and data
        """
        return self._do("GET", endpoint, params = params)
    



        



