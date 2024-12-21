#! /usr/bin/env python3
# Copyright (C) 2019-2025 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import json
import re
from collections import OrderedDict

class CodeMetaManipulator(object):
    def load(self, filename='codemeta.json'):
        with open(filename, 'rb') as fp:
            self.data = json.load(fp, object_pairs_hook=OrderedDict)

    @staticmethod
    def _dict_entry_cmp(dict1, dict2, field):
        if (field in dict1) and (field in dict2):
            return dict1[field] == dict2[field]
        else:
            return False

    @classmethod
    def find_person_entry(cls, person_list, matchdict):
        for entry in person_list:
            if cls._dict_entry_cmp(entry, matchdict, 'email'):
                return entry
            if cls._dict_entry_cmp(entry, matchdict, 'familyName') \
                    and cls._dict_entry_cmp(entry, matchdict, 'givenName'):
                return entry
        return None

    @staticmethod
    def update_person_entry(entry, matchdict):
        if entry is None:
            entry = OrderedDict()
            entry['@type'] = 'Person'
        
        # Aquí se cambian los nombres de los campos
        for field in ('givenName', 'familyName', 'email', 'affiliation'):
            val = matchdict.get(field, None)
            if val is not None:
                entry[field] = val
        
        # Cambiar 'orcid' a 'identifier' y asignarlo
        if matchdict.get('orcid'):
            entry['identifier'] = matchdict['orcid']
        
        return entry

    @staticmethod
    def is_valid_email(email):
        """ Verifica si el correo electrónico tiene un formato válido. """
        email_regex = r"(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)"
        return bool(re.match(email_regex, email))

    def handle_person_list_file(self, filename, cm_field_name):
        with open(filename, 'r', encoding='utf8') as fp:
            findregex = re.compile(
                r'^(?P<familyName>[-\w\s]*[-\w]),\s*'
                r'(?P<givenName>[-\w\s]*[-\w])\s*'
                r'(?:\[(?P<email>[^\]]+)\])?\s*'
                r'(?:\[(?P<orcid>https?://orcid.org/[^\]]+)\])?\s*'
                r'(?:\[(?P<affiliation>[^\]]+)\])?$'
            )
            person_list = self.data.setdefault(cm_field_name, [])
            for line in fp:
                line = line.strip()
                m = findregex.match(line)
                if m is None:
                    raise RuntimeError(f"Could not analyze line: {line!r}")
                matchdict = m.groupdict()

                # Verificar si el email es válido
                if matchdict.get('email'):
                    if not self.is_valid_email(matchdict['email']):
                        matchdict['email'] = None  # Asignar None si el correo es inválido

                found_entry = self.find_person_entry(person_list, matchdict)
                entry = self.update_person_entry(found_entry, matchdict)
                if found_entry is None:
                    person_list.append(entry)

    def save(self, filename='codemeta.json'):
        with open(filename, 'w', encoding='utf8') as fp:
            json.dump(self.data, fp, indent=2)
            fp.write('\n')


def main():
    parser = argparse.ArgumentParser(description='Update codemeta.json')
    parser.add_argument('--set-version', dest='newversion')
    args = parser.parse_args()

    cm = CodeMetaManipulator()
    cm.load()
    if args.newversion is not None:
        cm.data['softwareVersion'] = args.newversion
    cm.handle_person_list_file('AUTHORS', 'author')
    cm.handle_person_list_file('CONTRIBUTORS', 'contributor')
    cm.save()


if __name__ == '__main__':
    main()
