from typing import Any

from modflow_devtools.misc import try_get_enum_value


class Filters:
    def title(component_name: tuple[str, str]) -> str:
        """
        The input component's unique title.
        """
        l, r = component_name
        return f"{l}-{r}idm"

    @staticmethod
    def value(v: Any) -> str:
        """
        Format a value to appear in the RHS of an assignment or argument-passing expression.
        """
        v = try_get_enum_value(v)
        v = f"'{v}'"
        if isinstance(v, bool):
            v = ".true." if v else ".false."
        return v