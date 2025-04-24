from abc import ABC, abstractmethod
from exceptions import ValidationError
from enums.enums import FilterError
class BaseCollection(ABC):

    @abstractmethod
    def _get_items(self):
        pass

    def filter(self, **kwargs):
        """
        Filter the items in the collection by the given keyword arguments.
        :param kwargs: A dictionary of keyword arguments to filter the items by
        :return: A new Collection object containing the filtered items
        """
        invalid_attrs = self._get_invalid_attrs(kwargs)
        if invalid_attrs:
            raise ValidationError(f"Invalid attribute {invalid_attrs} for {self.filter.__name__}")
        return self._filter_items(kwargs)

    def _get_invalid_attrs(self, kwargs):
        valid_attrs = self._get_items()[0].__dict__.keys()
        return [k for k in kwargs if k not in valid_attrs]

    def _filter_items(self, kwargs):
        filtered = [item for item in self._get_items()
                   if all(getattr(item, k) == v 
                         for k, v in kwargs.items())]
        return self.__class__(data=filtered)