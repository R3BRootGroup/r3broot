#! /usr/bin/env python3
# Copyright (C) 2019-2024 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import re
import json
from collections import OrderedDict
from argparse import ArgumentParser
import os

class Manipulator(object):
    def __str__(self):
        return self.__class__.__name__

    def load(self, filename=None):
        if filename is None:
            filename = self.default_filename
        try:
            with open(filename, 'r', encoding='utf8') as fp:
                self.data = json.load(fp, object_pairs_hook=OrderedDict)
        except FileNotFoundError as e:
            print(f'Error: The file {filename} was not found. {e}')
        except json.JSONDecodeError as e:
            print(f'Error in JSON file: {e}')
            raise  # Raise the error and stop execution

    def save(self, filename=None, indent=2):
        if filename is None:
            filename = self.default_filename
        try:
            with open(filename, 'w', encoding='utf8') as fp:
                json.dump(self.data, fp, indent=indent, ensure_ascii=False)
                fp.write('\n')
        except IOError as e:
            print(f'Error saving the file {filename}: {e}')

    @staticmethod
    def _dict_entry_cmp(dict1, dict2, field1, field2=None):
        if field2 is None:
            field2 = field1
        if (field1 in dict1) and (field2 in dict2):
            return dict1[field1] == dict2[field2]
        return False

    def _handle_person_list_file(self, filename, field_name, **kwargs):
        # Check if the file exists before trying to open it
        if not os.path.exists(filename):
            print(f'File {filename} not found, skipping.')
            return

        try:
            with open(filename, 'r', encoding='utf8') as fp:
                person_list = self.data.setdefault(field_name, [])
                for i, line in enumerate(fp, start=0):
                    line = line.strip()
                    m = self.findregex.match(line)
                    if m is None:
                        print(f"Could not analyze line: {line}")
                        continue  # Skip malformed lines
                    found_entry = self._find_person_entry(person_list, m.groupdict())
                    entry = self.update_person_entry(found_entry, m.groupdict(), **kwargs)
                    if found_entry is None:
                        person_list.insert(i, entry)
        except FileNotFoundError as e:
            print(f'Error: The file {filename} was not found. {e}')
        except RuntimeError as e:
            print(e)

class ZenodoManipulator(Manipulator):
    default_filename = '.zenodo.json'
    # Improved regular expression to handle different author formats
    findregex = re.compile(r'^(?P<name>[-\w\s,]*[-\w])\s*'  # Name
                           r'(?:\[(?P<email>\S+@\S+)\])?\s*'  # Optional email
                           r'(\[https://orcid\.org/(?P<orcid>\S+)\])?\s*'  # Optional ORCID
                           r'(\[(?P<affiliation>[^\]]+)\])?$')  # Optional affiliation

    @classmethod
    def _find_person_entry(cls, person_list, matchdict):
        # First attempt to match using ORCID
        for entry in person_list:
            if cls._dict_entry_cmp(entry, matchdict, 'orcid'):
                return entry
        # If not, attempt to match by name
        for entry in person_list:
            if cls._dict_entry_cmp(entry, matchdict, 'name'):
                return entry
        return None

    @staticmethod
    def update_person_entry(entry, matchdict, contributor_type=None):
        if entry is None:
            entry = OrderedDict()
            if contributor_type:
                entry['type'] = contributor_type
        for field in ('name', 'orcid', 'email', 'affiliation'):
            val = matchdict.get(field, None)
            if val is not None:
                entry[field] = val
        return entry

    def update_authors(self):
        # Update the authors and contributors if the files exist
        self._handle_person_list_file('AUTHORS', 'creators')
        self._handle_person_list_file('CONTRIBUTORS', 'contributors', contributor_type='Other')

    def save(self, filename=None):
        super().save(filename, 4)

    def version(self, new_version):
        self.data['version'] = new_version


def main():
    parser = ArgumentParser(description='Update .zenodo.json')
    parser.add_argument('--set-version', dest='newversion')
    args = parser.parse_args()

    for manipulator in (ZenodoManipulator(),):
        try:
            manipulator.load()
        except Exception as e:
            print(f'*** Error loading {manipulator}: {e}')
            continue
        if args.newversion is not None:
            manipulator.version(args.newversion)
        manipulator.update_authors()
        manipulator.save()


if __name__ == '__main__':
    main()
