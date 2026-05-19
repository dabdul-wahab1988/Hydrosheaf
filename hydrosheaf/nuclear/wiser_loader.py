from __future__ import annotations

import math
import re
import zipfile
from pathlib import Path
from typing import Dict, Optional
from xml.etree.ElementTree import iterparse

import pandas as pd

from hydrosheaf.nuclear.input_history import InputHistory

class WiserInputLibrary:
    """Intelligent lookup for localized atmospheric tracer inputs from IAEA WISER."""

    REQUIRED_COLUMNS = (
        "Sample Date",
        "Sample Site Name",
        "Latitude",
        "Longitude",
        "WMO Code",
        "Measurand Symbol",
        "Measurand Amount",
        "Measurand Uncertainty",
    )
    XLSX_REQUIRED_COLS = {
        "L": "Sample Site Name",
        "M": "Latitude",
        "N": "Longitude",
        "P": "WMO Code",
        "T": "Sample Date",
        "AI": "Measurand Symbol",
        "AM": "Measurand Amount",
        "AN": "Measurand Uncertainty",
    }
    XLSX_ROW_RE = re.compile(br"<row\b[^>]*>.*?</row>")
    XLSX_CELL_RE = re.compile(br"<c\b(?![^>]*/>)(?P<attrs>[^>]*)>(?P<body>.*?)</c>")
    XLSX_REF_RE = re.compile(br'\br="([A-Z]+)\d+"')
    XLSX_VALUE_RE = re.compile(br"<v>(.*?)</v>")
    
    def __init__(self, data_path: str | Path):
        self.path = Path(data_path)
        if not self.path.exists():
            raise FileNotFoundError(f"WISER data not found at {self.path}")
        self._df: Optional[pd.DataFrame] = None
        self._shared_strings: Optional[list[str]] = None
        self._tracer_frame_cache: Dict[tuple[str, ...], pd.DataFrame] = {}
        self._history_cache: Dict[tuple[str, str], InputHistory] = {}
        self._station_frame_cache: Dict[str, pd.DataFrame] = {}

    @property
    def df(self) -> pd.DataFrame:
        """Load the configured WISER file lazily.

        The North America workbook is large enough that importing this module
        should not parse it. M3 explicitly uses the North America workbook only.
        """
        if self._df is None:
            self._df = self._load_frame()
        return self._df

    def _load_frame(self) -> pd.DataFrame:
        # Support both CSV and Excel
        if self.path.suffix.lower() == '.csv':
            df = pd.read_csv(self.path, encoding='latin-1', usecols=lambda c: c in self.REQUIRED_COLUMNS)
        elif self.path.suffix.lower() in ['.xlsx', '.xls']:
            df = self._read_xlsx_required_columns()
        else:
            raise ValueError(f"Unsupported file format: {self.path.suffix}")
            
        # Standardize date and year
        df['Sample Date'] = pd.to_datetime(df['Sample Date'], errors='coerce', utc=True)
        # Use simple mid-month year representation
        df['Year'] = df['Sample Date'].dt.year + (df['Sample Date'].dt.month - 0.5) / 12.0
        if "Latitude" in df.columns:
            df["Latitude"] = pd.to_numeric(df["Latitude"], errors="coerce")
        if "Longitude" in df.columns:
            df["Longitude"] = pd.to_numeric(df["Longitude"], errors="coerce")
        df["Measurand Amount"] = pd.to_numeric(df["Measurand Amount"], errors="coerce")
        df["Measurand Uncertainty"] = pd.to_numeric(df["Measurand Uncertainty"], errors="coerce")
        return df

    def _read_xlsx_required_columns(self) -> pd.DataFrame:
        """Read the WISER observations sheet without the heavy Excel engine.

        The workbook stores a single large XML sheet. Streaming the handful of
        required columns is much faster than asking pandas/openpyxl to build the
        full workbook object.
        """
        namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

        def col_number(cell_ref: str) -> int:
            letters = re.sub(r"\d+", "", cell_ref)
            number = 0
            for letter in letters:
                number = number * 26 + ord(letter.upper()) - 64
            return number

        def shared_strings(zf: zipfile.ZipFile) -> list[str]:
            strings: list[str] = []
            with zf.open("xl/sharedStrings.xml") as handle:
                for _, elem in iterparse(handle, events=("end",)):
                    if elem.tag == namespace + "si":
                        strings.append("".join(t.text or "" for t in elem.iter(namespace + "t")))
                        elem.clear()
            return strings

        def cell_value(cell, strings: list[str]):
            cell_type = cell.attrib.get("t")
            value = cell.find(namespace + "v")
            if cell_type == "s" and value is not None:
                return strings[int(value.text)]
            if cell_type == "inlineStr":
                return "".join(t.text or "" for t in cell.iter(namespace + "t"))
            return value.text if value is not None else None

        rows: list[dict[str, object]] = []
        with zipfile.ZipFile(self.path) as zf:
            strings = self._get_shared_strings(zf)
            with zf.open("xl/worksheets/sheet1.xml") as handle:
                target_cols: dict[int, str] = {}
                for _, elem in iterparse(handle, events=("end",)):
                    if elem.tag != namespace + "row":
                        continue
                    if not target_cols:
                        cells = {
                            col_number(cell.attrib["r"]): cell_value(cell, strings)
                            for cell in elem.findall(namespace + "c")
                        }
                        header_by_col = {idx: str(value) for idx, value in cells.items() if value is not None}
                        target_cols = {
                            idx: header
                            for idx, header in header_by_col.items()
                            if header in self.REQUIRED_COLUMNS
                        }
                    else:
                        row = {header: None for header in target_cols.values()}
                        for cell in elem.findall(namespace + "c"):
                            idx = col_number(cell.attrib["r"])
                            header = target_cols.get(idx)
                            if header is not None:
                                row[header] = cell_value(cell, strings)
                        if any(value is not None for value in row.values()):
                            rows.append(row)
                    elem.clear()
        return pd.DataFrame(rows, columns=list(self.REQUIRED_COLUMNS))

    def _get_shared_strings(self, zf: Optional[zipfile.ZipFile] = None) -> list[str]:
        """Return cached shared strings for the configured WISER workbook."""
        if self._shared_strings is not None:
            return self._shared_strings

        namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
        strings: list[str] = []
        close_zip = zf is None
        workbook = zf if zf is not None else zipfile.ZipFile(self.path)
        try:
            with workbook.open("xl/sharedStrings.xml") as handle:
                try:
                    from lxml import etree

                    iterator = etree.iterparse(handle, events=("end",), tag=namespace + "si")
                    for _, elem in iterator:
                        strings.append("".join(t.text or "" for t in elem.iter(namespace + "t")))
                        elem.clear()
                except ImportError:
                    for _, elem in iterparse(handle, events=("end",)):
                        if elem.tag == namespace + "si":
                            strings.append("".join(t.text or "" for t in elem.iter(namespace + "t")))
                            elem.clear()
        finally:
            if close_zip:
                workbook.close()

        self._shared_strings = strings
        return strings

    @staticmethod
    def _cell_text_from_row(row_xml: bytes, col: str, strings: list[str]) -> Optional[str]:
        col_bytes = col.encode("ascii")
        pattern = re.compile(
            br'<c\b(?![^>]*/>)[^>]*\br="' + col_bytes + br'\d+"[^>]*>.*?</c>'
        )
        match = pattern.search(row_xml)
        if match is None:
            return None
        cell_xml = match.group(0)
        value = re.search(br"<v>(.*?)</v>", cell_xml)
        if value is None:
            return None
        raw = value.group(1).decode("utf-8")
        tag = cell_xml.split(b">", 1)[0]
        if b't="s"' in tag:
            return strings[int(raw)]
        return raw

    @classmethod
    def _cell_texts_from_row(cls, row_xml: bytes, strings: list[str]) -> dict[str, Optional[str]]:
        values: dict[str, Optional[str]] = {header: None for header in cls.XLSX_REQUIRED_COLS.values()}
        for match in cls.XLSX_CELL_RE.finditer(row_xml):
            attrs = match.group("attrs")
            ref = cls.XLSX_REF_RE.search(attrs)
            if ref is None:
                continue
            col = ref.group(1).decode("ascii")
            header = cls.XLSX_REQUIRED_COLS.get(col)
            if header is None:
                continue
            value = cls.XLSX_VALUE_RE.search(match.group("body"))
            if value is None:
                continue
            raw = value.group(1).decode("utf-8")
            values[header] = strings[int(raw)] if b't="s"' in attrs else raw
        return values

    @staticmethod
    def _shared_cell_re(col: str, ids: set[int]) -> Optional[re.Pattern[bytes]]:
        if not ids:
            return None
        col_bytes = col.encode("ascii")
        id_alt = b"|".join(str(idx).encode("ascii") for idx in sorted(ids))
        return re.compile(
            br'<c\b[^>]*\br="' + col_bytes + br'\d+"[^>]*>\s*<v>(?:'
            + id_alt
            + br")</v>\s*</c>"
        )

    def _read_xlsx_filtered_rows(self, station_name_or_wmo: str, symbols: list[str]) -> pd.DataFrame:
        """Fast path for station/tracer lookup in the large WISER workbook.

        The North America WISER file is a single large worksheet. For M3 we
        repeatedly need one localized input history, not the full observation
        table, so this scans row XML and materializes only matching rows.
        """
        station_key = str(station_name_or_wmo).upper()
        symbol_keys = {str(symbol).upper() for symbol in symbols}
        rows: list[dict[str, object]] = []

        with zipfile.ZipFile(self.path) as zf:
            strings = self._get_shared_strings(zf)
            station_ids = {
                idx
                for idx, text in enumerate(strings)
                if station_key in str(text).upper()
            }
            symbol_ids = {
                idx
                for idx, text in enumerate(strings)
                if str(text).upper() in symbol_keys
            }
            if not symbol_ids:
                return pd.DataFrame(columns=list(self.REQUIRED_COLUMNS))

            station_cell = self._shared_cell_re("L", station_ids)
            symbol_cell = self._shared_cell_re("AI", symbol_ids)
            sheet_xml = zf.read("xl/worksheets/sheet1.xml")

        if station_cell is not None:
            for station_match in station_cell.finditer(sheet_xml):
                row_start = sheet_xml.rfind(b"<row", 0, station_match.start())
                row_end = sheet_xml.find(b"</row>", station_match.end())
                if row_start < 0 or row_end < 0:
                    continue
                row_xml = sheet_xml[row_start:row_end + len(b"</row>")]
                if symbol_cell is not None and symbol_cell.search(row_xml) is not None:
                    rows.append(self._cell_texts_from_row(row_xml, strings))
        else:
            for row_match in self.XLSX_ROW_RE.finditer(sheet_xml):
                row_xml = row_match.group(0)
                if symbol_cell is None or symbol_cell.search(row_xml) is None:
                    continue

                station_value = self._cell_text_from_row(row_xml, "L", strings)
                wmo_value = self._cell_text_from_row(row_xml, "P", strings)
                station_ok = station_key in str(station_value or "").upper() or station_key in str(wmo_value or "").upper()
                if station_ok:
                    rows.append(self._cell_texts_from_row(row_xml, strings))

        return pd.DataFrame(rows, columns=list(self.REQUIRED_COLUMNS))

    def _prepare_observation_frame(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["Sample Date"] = pd.to_datetime(df["Sample Date"], errors="coerce", utc=True)
        df["Year"] = df["Sample Date"].dt.year + (df["Sample Date"].dt.month - 0.5) / 12.0
        for col in ("Latitude", "Longitude", "Measurand Amount", "Measurand Uncertainty"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df

    def _read_xlsx_symbol_rows(self, symbols: list[str]) -> pd.DataFrame:
        symbol_keys = {str(symbol).upper() for symbol in symbols}
        cache_key = tuple(sorted(symbol_keys))
        cached = self._tracer_frame_cache.get(cache_key)
        if cached is not None:
            return cached

        rows: list[dict[str, object]] = []
        with zipfile.ZipFile(self.path) as zf:
            strings = self._get_shared_strings(zf)
            symbol_ids = {
                idx
                for idx, text in enumerate(strings)
                if str(text).upper() in symbol_keys
            }
            if not symbol_ids:
                empty = pd.DataFrame(columns=list(self.REQUIRED_COLUMNS))
                self._tracer_frame_cache[cache_key] = empty
                return empty
            symbol_cell = self._shared_cell_re("AI", symbol_ids)
            sheet_xml = zf.read("xl/worksheets/sheet1.xml")

        if symbol_cell is not None:
            for symbol_match in symbol_cell.finditer(sheet_xml):
                row_start = sheet_xml.rfind(b"<row", 0, symbol_match.start())
                row_end = sheet_xml.find(b"</row>", symbol_match.end())
                if row_start < 0 or row_end < 0:
                    continue
                row_xml = sheet_xml[row_start:row_end + len(b"</row>")]
                rows.append(self._cell_texts_from_row(row_xml, strings))

        frame = self._prepare_observation_frame(pd.DataFrame(rows, columns=list(self.REQUIRED_COLUMNS)))
        self._tracer_frame_cache[cache_key] = frame
        return frame

    def _history_from_subset(self, subset: pd.DataFrame, measurand: str) -> InputHistory:
        if subset.empty:
            raise ValueError(f"No data found for {measurand} in {self.path.name}")

        subset = subset.dropna(subset=['Year', 'Measurand Amount'])
        if subset.empty:
            raise ValueError(f"No finite data found for {measurand} in {self.path.name}")
        subset = subset.sort_values('Year')
        daily = subset.groupby('Year').agg({
            'Measurand Amount': 'mean',
            'Measurand Uncertainty': 'mean'
        }).reset_index()

        if measurand.lower() == "tritium":
            pre_bomb = pd.DataFrame({
                'Year': [1900.0, 1910.0, 1920.0, 1930.0, 1940.0, 1950.0],
                'Measurand Amount': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
                'Measurand Uncertainty': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
            })
            daily = pd.concat([pre_bomb, daily]).sort_values('Year').reset_index(drop=True)

        return InputHistory(
            daily['Year'].values,
            daily['Measurand Amount'].values,
            daily['Measurand Uncertainty'].fillna(0.0).values
        )

    @staticmethod
    def _symbols_for_measurand(measurand: str) -> list[str]:
        if measurand.lower() == "tritium":
            return ["Tritium", "H3", "3H"]
        return [measurand]

    def get_nearest_input_history(self, latitude: float, longitude: float, measurand: str = "Tritium") -> InputHistory:
        """Retrieve the nearest WISER station history for a sample coordinate."""
        lat = float(latitude)
        lon = float(longitude)
        if not pd.notna(lat) or not pd.notna(lon):
            raise ValueError("Latitude and longitude are required for nearest WISER lookup.")

        cache_key = (f"nearest:{round(lat, 3)}:{round(lon, 3)}", str(measurand).lower())
        cached = self._history_cache.get(cache_key)
        if cached is not None:
            return cached

        symbols = self._symbols_for_measurand(measurand)
        if self.path.suffix.lower() == ".xlsx" and self._df is None:
            tracer_frame = self._read_xlsx_symbol_rows(symbols)
        else:
            tracer_frame = self.df[self.df["Measurand Symbol"].isin(symbols)]
        station_frame_key = str(measurand).lower()
        stations = self._station_frame_cache.get(station_frame_key)
        if stations is None:
            stations = (
                tracer_frame
                .dropna(subset=["Latitude", "Longitude"])
                .groupby(["Sample Site Name", "WMO Code"], dropna=False)
                .agg(Latitude=("Latitude", "mean"), Longitude=("Longitude", "mean"), n=("Measurand Amount", "count"))
                .reset_index()
            )
            stations = stations[stations["n"] > 0]
            self._station_frame_cache[station_frame_key] = stations
        if stations.empty:
            raise ValueError(f"No coordinate-indexed WISER stations found for {measurand} in {self.path.name}")

        lat_scale = math.cos(math.radians(lat))
        dlat = stations["Latitude"] - lat
        dlon = (stations["Longitude"] - lon) * lat_scale
        nearest = stations.loc[(dlat * dlat + dlon * dlon).idxmin()]
        station_cache_key = (
            f"station:{nearest['Sample Site Name']}:{nearest['WMO Code']}",
            str(measurand).lower(),
        )
        station_cached = self._history_cache.get(station_cache_key)
        if station_cached is not None:
            self._history_cache[cache_key] = station_cached
            return station_cached

        station_mask = pd.Series(False, index=tracer_frame.index)
        station_name = nearest["Sample Site Name"]
        wmo_code = nearest["WMO Code"]
        if pd.notna(station_name):
            station_mask |= tracer_frame["Sample Site Name"] == station_name
        if pd.notna(wmo_code):
            station_mask |= tracer_frame["WMO Code"].astype(str) == str(wmo_code)
        subset = tracer_frame[station_mask]
        history = self._history_from_subset(subset, measurand)
        self._history_cache[station_cache_key] = history
        self._history_cache[cache_key] = history
        return history
        
    def get_site_chemistry(self, station_name_or_wmo: str) -> Dict[str, float]:
        """Retrieve average geochemical parameters for a site to inform A0 corrections."""
        # Standardize site name search
        subset = self.df[
            (self.df['Sample Site Name'].str.contains(station_name_or_wmo, case=False, na=False)) | 
            (self.df['WMO Code'].astype(str).str.contains(str(station_name_or_wmo), na=False))
        ]
        
        if subset.empty:
            return {}

        # Look for the keys we need
        chem_map = {
            'C13': 'd13c_permil',
            'HCO3': 'HCO3_mmoll',
            'pH': 'pH'
        }
        
        result = {}
        for symbol, target_key in chem_map.items():
            vals = subset[subset['Measurand Symbol'] == symbol]['Measurand Amount']
            if not vals.empty:
                result[target_key] = float(vals.mean())
        
        return result

    def get_input_history(self, station_name_or_wmo: str, measurand: str = "Tritium") -> InputHistory:
        """Retrieve history for a specific station and tracer."""
        cache_key = (str(station_name_or_wmo).lower(), str(measurand).lower())
        cached = self._history_cache.get(cache_key)
        if cached is not None:
            return cached

        # Symbol mapping (IAEA WISER uses different codes)
        symbols = self._symbols_for_measurand(measurand)

        if self.path.suffix.lower() == ".xlsx" and self._df is None:
            subset = self._read_xlsx_filtered_rows(station_name_or_wmo, symbols)
            subset = self._prepare_observation_frame(subset)
        else:
            # Filter for station
            subset = self.df[
                (self.df['Sample Site Name'].str.contains(station_name_or_wmo, case=False, na=False)) | 
                (self.df['WMO Code'].astype(str).str.contains(str(station_name_or_wmo), na=False))
            ]
        
            # Filter for tracer
            subset = subset[subset['Measurand Symbol'].isin(symbols)]
        
        if subset.empty:
            # Try searching by Country if station not found? 
            # For now just raise to be strict.
            raise ValueError(f"No data found for {station_name_or_wmo} / {measurand} in {self.path.name}")
            
        history = self._history_from_subset(subset, measurand)
        self._history_cache[cache_key] = history
        return history

# Regional instances for intelligent automated selection
# For North American USGS validation
NA_WISER_PATH = Path(__file__).resolve().parents[2] / "data" / "wiser_north_america.xlsx"
WISER_NA = WiserInputLibrary(NA_WISER_PATH) if NA_WISER_PATH.exists() else None

# M3 and USGS validation must not implicitly fall back to the full/global WISER
# data file. Keep this symbol for backwards-compatible imports, but do not
# initialize it.
GLOBAL_WISER_PATH = None
WISER_GLOBAL = None

# Default legacy pointer
WISER_LIB = WISER_NA
